# -*- coding: utf-8 -*-
"""
Core pipeline functions for HobIT: housekeeping, the TBLASTN search, parsing
the XML into Hit objects, clustering them into MergedHits, and enriching those
with genome windows and query sequences.

Here define:
    serializer(args, obj, filename)
    unserializer(args, pickle_path)
    progress(args, pct)
    file_sys_prep(args)
    blast_cmdline_checker()
    _blast_db_exists(genome_path)
    blast_executer(args)
    get_blast_hits(args)
    _norm_coords(hit)
    _hit_sort_key(hit)
    make_mergedhits(args, all_hits)
    __get_locations_in_genome(mergedhits)
    __retrieve_genome_chunks(args, locations, genome_dict)
    __attach_genome_chunks(mergedhits, genome_chunks)
    __get_best_hitting_queries(mergedhits)
    __retrieve_query_sequences(args, query_ids)
    __attach_query_sequences(mergedhits, query_sequences)
    enrich_mergedhits(args, mergedhits, genome_dict)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
Adapted from the unpublished HobIT prototype by Thomas Barton-Owen.
"""

import pickle
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIXML

from hobit.bioinformatics_tools import BioinformaticsTools
from hobit.structures import Hit, MergedHit


# FUNCS FOR HOUSEKEEPING

def serializer(args, obj, filename):
    # pickle the object into the output folder, so a run can be resumed with -pic
    filename = Path(args.out).joinpath(filename)

    with open(filename, "wb") as outfile:
        pickle.dump(obj, outfile)

    if args.verbose:
        print(f"Saved to {filename}")


def unserializer(args, pickle_path):
    # load a previously pickled pipeline state
    with open(pickle_path, "rb") as infile:
        return pickle.load(infile)


def progress(args, pct):
    if args.verbose:
        print(f"Progress: {pct}%")


def file_sys_prep(args):
    # check the input files exist and create the output folder; quit if not
    quit_bool = False
    genome = args.genome if isinstance(args.genome, str) else args.genome[0]

    if not Path(genome).is_file():
        print(f"Genome file not found: {genome}")
        quit_bool = True

    if not Path(args.query).is_file():
        print(f"Query file not found: {args.query}")
        quit_bool = True

    if quit_bool:
        sys.exit(1)

    Path(args.out).mkdir(parents=True, exist_ok=True)


# FUNCS FOR HANDLING BLAST SEARCHES

def blast_cmdline_checker():
    # returns True if tblastn is on the PATH
    try:
        return subprocess.run(["tblastn", "-version"],
                              capture_output=True, text=True).returncode == 0
    except FileNotFoundError:
        return False


def _blast_db_exists(genome_path):
    # makeblastdb writes its index next to the FASTA with the extension appended
    # (genome.fa -> genome.fa.nhr); .nhr is the legacy index, .nsi/.nsq are v5
    genome_path = Path(genome_path)
    #
    for suffix in (".nhr", ".nsi", ".nsq"):
        index_path = genome_path.with_name(genome_path.name + suffix)
        if index_path.exists():
            return True

    return False


def blast_executer(args):
    # run TBLASTN (protein queries vs nucleotide genome), building the database if needed
    xml_path = Path(args.out) / "blast_output.xml"

    if args.verbose:
        print("Running TBLASTN search...")

    if not blast_cmdline_checker():
        sys.exit("TBLASTN not found in PATH. Please install BLAST+.")

    genome = str(args.genome) if isinstance(args.genome, (str, Path)) else str(args.genome[0])

    # checking for an existing database stops makeblastdb re-running every time
    if not args.no_mkblastdb and not _blast_db_exists(genome):
        if args.verbose:
            print(f"Creating BLAST database for {genome}...")
        subprocess.call(["makeblastdb", "-in", genome,
                         "-dbtype", "nucl", "-parse_seqids"])

    subprocess.call([
        "tblastn",
        "-db", genome,
        "-query", args.query,
        "-outfmt", "5",                  # 5 = XML, which is what NCBIXML reads
        "-out", str(xml_path),
        "-num_threads", str(args.threads)
    ])

    args.xml = str(xml_path)

    return args


# FUNCS FOR READING THE XML INTO HIT OBJECTS

def get_blast_hits(args):
    if args.verbose:
        print("Processing BLAST results...")

    with open(args.xml) as handle:
        blast_records = list(NCBIXML.parse(handle))

    all_hits = []

    for record in blast_records:
        # names are sanitised here so they match the keys SeqIO will produce later
        query_name = BioinformaticsTools.sanitise(record.query)

        for alignment in record.alignments:
            scaffold = BioinformaticsTools.sanitise(alignment.hit_id or alignment.hit_def)

            for hsp in alignment.hsps:
                # subject coordinates are stored exactly as BLAST gave them,
                # so on the minus strand sbjct_start > sbjct_end
                all_hits.append(Hit(
                    query=query_name,
                    query_location=(hsp.query_start, hsp.query_end),
                    score=hsp.score,
                    evalue=hsp.expect,
                    align_list=[hsp.query, hsp.match, hsp.sbjct],
                    scaffold=scaffold,
                    scaffold_location=(hsp.sbjct_start, hsp.sbjct_end),
                    frame=hsp.frame[1] if hsp.frame else 1   # the sign carries the strand
                ))

    if args.verbose:
        print(f"    {len(all_hits)} discrete hits")

    return all_hits


# FUNCS FOR CLUSTERING HITS INTO MERGEDHITS

def _norm_coords(hit):
    # (leftmost, rightmost) regardless of which way BLAST wrote them
    start, end = hit.scaffold_location

    return min(start, end), max(start, end)


def _hit_sort_key(hit):
    # sort by scaffold, then strand, then leftmost coordinate
    strand = "+" if hit.frame > 0 else "-"
    leftmost, _ = _norm_coords(hit)

    return (hit.scaffold, strand, leftmost)


def make_mergedhits(args, all_hits):
    if args.verbose:
        print(f"Merging hits... [cluster = {args.cluster} nucleotides]")

    if not all_hits:
        return []

    # after this sort, members of one cluster always sit next to each other
    all_hits.sort(key=_hit_sort_key)

    merged_hits = []
    current_hits = [all_hits[0]]
    _, cluster_end = _norm_coords(all_hits[0])

    for next_hit in all_hits[1:]:
        # the cluster's first hit carries the scaffold/strand for comparison
        first_hit = current_hits[0]
        same_strand = (first_hit.frame > 0) == (next_hit.frame > 0)

        # distance is measured from the furthest point the cluster has reached,
        # not from the previous hit, so a long hit cannot make the cluster creep
        next_start, next_end = _norm_coords(next_hit)
        close_enough = (next_start - cluster_end) < args.cluster

        # requiring the same strand stops two genes facing opposite ways from merging
        if next_hit.scaffold == first_hit.scaffold and same_strand and close_enough:
            current_hits.append(next_hit)
            cluster_end = max(cluster_end, next_end)
        else:
            # too far away or a different scaffold/strand: close this cluster
            merged_hits.append(MergedHit(current_hits))
            current_hits = [next_hit]
            cluster_end = next_end

    merged_hits.append(MergedHit(current_hits))

    if args.verbose:
        print(f"    {len(merged_hits)} merged hits")

    return merged_hits


# FUNCS FOR ENRICHING THE MERGEDHITS WITH GENOMIC DATA

def __get_locations_in_genome(mergedhits):
    # build {scaffold: [leftmost position, ...]}
    # min() gives the leftmost coordinate on both strands, used as the window centre
    locations = {}

    for mh in mergedhits:
        locations.setdefault(mh.scaffold, []).append(min(mh.scaffold_location_guess))

    return locations


def __retrieve_genome_chunks(args, locations, genome_dict):
    # cut a +-args.window window from the FORWARD strand around each location
    # result: {scaffold: {location: (sequence, offset)}}, offset = 0-based window start
    genome_chunks = {scaffold: {} for scaffold in locations}

    for scaffold, scaffold_locations in locations.items():
        record = genome_dict.get(scaffold)

        if not record:
            if args.verbose:
                print(f"  Warning: scaffold {scaffold} not found in genome")
            continue

        for loc in scaffold_locations:
            start = max(0, loc - args.window)
            end = min(len(record.seq), loc + args.window)

            if args.verbose:
                # printed 1-based, because that is what a genome browser shows
                print(f"  {scaffold}:{loc} (nt) -> window {start + 1}-{end}")

            genome_chunks[scaffold][loc] = (str(record.seq[start:end]), start)

    return genome_chunks


def __attach_genome_chunks(mergedhits, genome_chunks):
    # hand each MergedHit its genome window, keyed on the same leftmost coordinate
    for mh in mergedhits:
        chunk = genome_chunks.get(mh.scaffold, {}).get(min(mh.scaffold_location_guess))

        if chunk:
            mh.attach_genome_sequence(chunk[0], chunk[1])
        else:
            # no window found: leave the hit intact but empty, so the pipeline
            # carries on and the failure shows up as a skipped annotation
            mh.attach_genome_sequence("", 0)

    return mergedhits


# FUNCS FOR ENRICHING THE MERGEDHITS WITH QUERY DATA

def __get_best_hitting_queries(mergedhits):
    # set of query IDs that produced a best hit somewhere
    return {mh.best_query for mh in mergedhits}


def __retrieve_query_sequences(args, query_ids):
    # read the protein sequences of those queries out of the query file
    sequences = {}

    for record in SeqIO.parse(args.query, "fasta"):
        query_id = BioinformaticsTools.sanitise(record.id)
        if query_id in query_ids:
            sequences[query_id] = str(record.seq)

    return sequences


def __attach_query_sequences(mergedhits, query_sequences):
    # hand each MergedHit the full protein sequence of its best query
    for mh in mergedhits:
        mh.attach_best_query_sequence(query_sequences.get(mh.best_query, ""))

    return mergedhits


# MASTER CONTROL FUNCTION FOR ENRICHING THE MERGEDHIT OBJECTS

def enrich_mergedhits(args, mergedhits, genome_dict):
    # sequence retrieval only: LPORF, HMM refinement and Exonerate are driven
    # from hobit_main, because their order matters and some are optional
    locations = __get_locations_in_genome(mergedhits)

    if args.verbose:
        print(f"Retrieving genome sequences... [window = +-{args.window} nucleotides]")

    genome_chunks = __retrieve_genome_chunks(args, locations, genome_dict)
    mergedhits = __attach_genome_chunks(mergedhits, genome_chunks)

    query_ids = __get_best_hitting_queries(mergedhits)
    query_sequences = __retrieve_query_sequences(args, query_ids)
    mergedhits = __attach_query_sequences(mergedhits, query_sequences)

    return mergedhits
