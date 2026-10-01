# -*- coding: utf-8 -*-
"""
Gene structure prediction with Exonerate. Every MergedHit is sent to Exonerate
in parallel; each locus gets up to three attempts anchored on its best BLAST
hits in turn, and each attempt widens the genome window until a result
validates.

Here define:
    class ExonerateResult
        __init__()
        get_full_gene_sequence()
        get_exon_count()
        get_intron_count()
    check_exonerate_installed()
    run_exonerate_annotation(args, mergedhits, genome_records=None)
    _run_exonerate_with_retries(args, mergedhit, genome_records, max_candidates=3)
    _run_exonerate_for_hit(args, mergedhit, genome_records)
    _parse_exonerate_output(output_text, genome_sequence)
    _infer_introns(exons, genome_sequence)
    _preliminary_validate(mergedhit, exonerate_result)
    _validate_exonerate_result(mergedhit, exonerate_result, args)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
"""

import difflib
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO

from hobit.bioinformatics_tools import BioinformaticsTools

# helix 3 of the homeodomain; more than one copy in a prediction means
# Exonerate has run two neighbouring genes together
HELIX3_PATTERN = re.compile(r"WF[QN]NRR", re.IGNORECASE)

# window sizes tried in turn when the first attempt finds nothing
WINDOW_LADDER = (2000, 5000, 10000)


# THE EXONERATE RESULT CONTAINER

class ExonerateResult:

    def __init__(self):
        # coordinates are 1-based, forward strand and relative to the genome window,
        # until gff_outputter adds the window offset
        self.exons = []               # (start, end, seq, phase) tuples, genomic order
        self.introns = []             # (start, end, seq, phase) tuples, genomic order
        self.cds_sequence = ""        # translational order (already reverse complemented on -)
        self.protein_sequence = ""
        self.score = 0
        self.gene_start = None
        self.gene_end = None
        self.strand = "+"

    def get_full_gene_sequence(self):
        # exons and introns put back together in genomic order
        if not self.exons:
            return ""

        features = [(start, seq) for start, _, seq, _ in self.exons]
        features += [(start, seq) for start, _, seq, _ in self.introns]
        features.sort(key=lambda f: f[0])

        return "".join(seq for _, seq in features)

    def get_exon_count(self):
        return len(self.exons)

    def get_intron_count(self):
        return len(self.introns)


# INSTALLATION CHECK

def check_exonerate_installed():
    path = shutil.which("exonerate")
    if not path:
        return False

    try:
        result = subprocess.run([path, "--version"], capture_output=True, text=True)
    except OSError:
        return False

    # some builds print the banner to stderr, hence checking both streams
    if "exonerate" in (result.stdout + result.stderr).lower():
        return True

    return result.returncode == 0


# MAIN ENTRY POINT

def run_exonerate_annotation(args, mergedhits, genome_records=None):
    if not check_exonerate_installed():
        if args.verbose:
            print("Exonerate not found in PATH. Skipping gene structure prediction.")
        return mergedhits

    if args.verbose:
        print("Predicting gene structures with Exonerate...")
        print(f"    Model: {args.exonerate_model}")

    # full scaffold records are needed to widen the window past the pre-extracted one
    if genome_records is None:
        genome_records = {}
        if any(mh.genome_sequence for mh in mergedhits):
            genome = str(args.genome) if isinstance(args.genome, (str, Path)) else str(args.genome[0])
            genome_records = {BioinformaticsTools.sanitise(rec.id): rec
                              for rec in SeqIO.parse(genome, "fasta")}
            if args.verbose:
                print(f"    Loaded {len(genome_records)} scaffolds for window expansion.")

    annotated_count = 0
    rejected_count = 0
    failed_count = 0

    max_workers = max(1, min(args.threads, len(mergedhits)))
    if max_workers > 1 and args.verbose:
        print(f"    Using {max_workers} parallel threads...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_hit = {}

        for mh in mergedhits:
            # a hit with no window or no query is counted as failed rather than raising
            if not mh.genome_sequence or not mh.best_query_sequence:
                failed_count += 1
                continue
            future = executor.submit(_run_exonerate_with_retries, args, mh, genome_records)
            future_to_hit[future] = mh

        completed = 0
        for future in as_completed(future_to_hit):
            mh = future_to_hit[future]
            completed += 1

            if args.verbose and completed % 10 == 0:
                print(f"    Progress: {completed}/{len(future_to_hit)}")

            # one locus blowing up must not abort the other several thousand
            try:
                result, status = future.result()
            except Exception:
                mh.notes.append("EXONERATE_FAILED")
                failed_count += 1
                continue

            if status == "ok":
                mh.exonerate_result = result
                mh.exon_count = result.get_exon_count()
                mh.intron_count = result.get_intron_count()
                mh.cds_sequence = result.cds_sequence
                mh.full_gene_sequence = result.get_full_gene_sequence()

                # note note note
                note = f"EXONERATE: {mh.exon_count} exon"
                if mh.exon_count != 1:
                    note += "s"
                if mh.intron_count > 0:
                    note += f", {mh.intron_count} intron"
                    if mh.intron_count != 1:
                        note += "s"
                mh.notes.append(note + f" (Score: {result.score})")
                annotated_count += 1

                # only replace improved_sequence when the prediction is actually better:
                # fewer stop codons always wins, on a tie the longer sequence does
                if result.protein_sequence:
                    current = mh.improved_sequence or ""
                    current_stops = current.count("*")
                    new_stops = result.protein_sequence.count("*")
                    current_len = len(current.replace("*", "").replace("-", ""))
                    new_len = len(result.protein_sequence.replace("*", "").replace("-", ""))

                    fewer_stops = new_stops < current_stops
                    same_stops_but_longer = new_stops == current_stops and new_len > current_len

                    if fewer_stops or same_stops_but_longer:
                        mh.improved_sequence = result.protein_sequence

            elif status == "rejected":
                mh.notes.append("EXONERATE_REJECTED")
                rejected_count += 1
            else:
                failed_count += 1

    if args.verbose:
        print(f"    Successfully annotated {annotated_count}/{len(mergedhits)} merged hits")
        if rejected_count:
            print(f"    Rejected: {rejected_count} (validation failed)")
        if failed_count:
            print(f"    Failed: {failed_count}")

    return mergedhits


# RETRY DISPATCHER

def _run_exonerate_with_retries(args, mergedhit, genome_records, max_candidates=3):
    # try the locus with each of its top BLAST hits until one validates;
    # returns (result or None, "ok" | "rejected" | "failed")

    # the MergedHit is temporarily repointed at each candidate, so remember its state
    original_state = {
        "scaffold": mergedhit.scaffold,
        "scaffold_location_guess": mergedhit.scaffold_location_guess,
        "frame": mergedhit.frame,
        "best_align_list": mergedhit.best_align_list,
        "best_align_query_location": mergedhit.best_align_query_location,
        "best_evalue": mergedhit.best_evalue,
        "best_score": mergedhit.best_score,
        "sequence": mergedhit.sequence,
        "improved_sequence": mergedhit.improved_sequence,
        "best_query": mergedhit.best_query,
    }

    candidates = sorted(mergedhit.hits, key=lambda h: h.evalue)[:max_candidates]
    any_rejected = False

    for rank, hit in enumerate(candidates, start=1):
        mergedhit.scaffold = hit.scaffold
        mergedhit.scaffold_location_guess = hit.scaffold_location
        mergedhit.frame = hit.frame
        mergedhit.best_align_list = hit.align_list
        mergedhit.best_align_query_location = hit.query_location
        mergedhit.best_evalue = hit.evalue
        mergedhit.best_score = hit.score
        mergedhit.sequence = hit.align_list[2].replace("-", "")
        mergedhit.improved_sequence = mergedhit.sequence
        mergedhit.best_query = hit.query

        if args.verbose:
            print(f"    [Candidate #{rank}] {mergedhit.scaffold}:"
                  f"{min(mergedhit.scaffold_location_guess)} (e={hit.evalue:.2e})")

        result = _run_exonerate_for_hit(args, mergedhit, genome_records)

        if result and result.exons:
            if _validate_exonerate_result(mergedhit, result, args):
                if rank > 1:
                    mergedhit.notes.append(f"EXONERATE_SECONDARY_HIT_USED (hit_rank={rank})")
                return result, "ok"

            # Exonerate produced something, but validation threw it out
            any_rejected = True

    # every candidate failed: restore the original state so no trace is left
    for key, value in original_state.items():
        setattr(mergedhit, key, value)

    return None, "rejected" if any_rejected else "failed"


# SINGLE HIT RUNNER

def _run_exonerate_for_hit(args, mergedhit, genome_records):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as query_file:
        query_file.write(f">query\n{mergedhit.best_query_sequence}\n")
        query_path = query_file.name

    # start from the window extracted during enrichment, then try each larger size;
    # a gene whose introns run past the window edge only appears once it is covered
    base_window = args.window
    windows = [base_window] + [w for w in WINDOW_LADDER if w > base_window]
    result = None

    try:
        for window in windows:
            if window == base_window:
                sequence = mergedhit.genome_sequence
                offset = mergedhit.genome_window_offset
            else:
                record = genome_records.get(mergedhit.scaffold)
                if not record:
                    break

                # min() is the leftmost position on either strand, so the wider window stays centred
                loc = min(mergedhit.scaffold_location_guess)
                offset = max(0, loc - window)
                window_end = min(len(record.seq), loc + window)
                sequence = str(record.seq[offset:window_end])

            if not sequence:
                break

            # the target is always the raw FORWARD-strand window: protein2genome searches
            # both strands itself, and a reverse complement would break the coordinates
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as target_file:
                target_file.write(f">target\n{sequence}\n")
                target_path = target_file.name

            try:
                exonerate = subprocess.run(
                    ["exonerate",
                     "--model", args.exonerate_model,
                     "--showtargetgff", "yes",
                     "--percent", "20",
                     "--bestn", "1",
                     query_path, target_path],
                    capture_output=True, text=True, timeout=300
                )

                stdout_lower = exonerate.stdout.lower()

                if "exon" in stdout_lower or "gene" in stdout_lower:
                    candidate = _parse_exonerate_output(exonerate.stdout, sequence)

                    if candidate.exons and _preliminary_validate(mergedhit, candidate):
                        result = candidate
                        # keep the window the coordinates were actually measured in
                        mergedhit.genome_sequence = sequence
                        mergedhit.genome_window_offset = offset
                        break

            except subprocess.TimeoutExpired:
                # a single slow locus should not stall the whole run
                if args.verbose:
                    print(f"    Warning: Exonerate timed out for {mergedhit.scaffold}")

            finally:
                Path(target_path).unlink(missing_ok=True)

    finally:
        Path(query_path).unlink(missing_ok=True)

    return result


# GFF OUTPUT PARSER

def _parse_exonerate_output(output_text, genome_sequence):
    # coordinates are kept as Exonerate gives them: 1-based, forward strand, window-relative
    result = ExonerateResult()

    gff_pattern = re.compile(
        r"(\S+)\s+exonerate:\S+\s+"
        r"(exon|intron|cds|gene|similarity|mrna)\s+"
        r"(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+([\.\d])\s+(.*)"
    )

    exons = []
    introns = []
    cds_segments = []

    for line in output_text.split("\n"):
        match = gff_pattern.search(line)
        if not match:
            continue

        _, feature, start_str, end_str, score_str, strand, phase, _ = match.groups()
        start, end = int(start_str), int(end_str)

        if feature == "gene":
            try:
                result.score = float(score_str)
            except ValueError:
                result.score = 0
            result.strand = strand
            continue

        # drop anything reaching outside the window; slicing it would silently
        # give a truncated sequence with correct-looking coordinates
        if start < 1 or end > len(genome_sequence):
            continue

        segment = (start, end, genome_sequence[start - 1:end], phase)

        if feature == "exon":
            exons.append(segment)
        elif feature == "cds":
            cds_segments.append(segment)
        elif feature == "intron":
            introns.append(segment)

    # where both cds and exon features exist, the cds boundaries are more reliable
    if cds_segments:
        exons = cds_segments

    exons.sort(key=lambda x: x[0])
    introns.sort(key=lambda x: x[0])

    result.exons = exons
    result.introns = introns

    if exons:
        result.gene_start = exons[0][0]
        result.gene_end = exons[-1][1]

        # concatenate the CDS segments in ascending genomic order
        raw_cds = "".join(seq for _, _, seq, _ in exons)

        # on the minus strand the whole CDS is reverse complemented, so the phase
        # comes from the LAST exon (the first one transcribed)
        if result.strand == "-":
            result.cds_sequence = BioinformaticsTools.reverse_complement(raw_cds)
            phase_str = exons[-1][3]
        else:
            result.cds_sequence = raw_cds
            phase_str = exons[0][3]

        # phase = how many bases to drop before the first whole codon; "." and "0" mean none
        phase_skip = int(phase_str) if phase_str.isdigit() else 0
        result.protein_sequence = BioinformaticsTools.translate(result.cds_sequence[phase_skip:], 1)

    # Exonerate does not always emit intron lines; derive them from the gaps
    if len(result.exons) > 1 and not result.introns:
        result.introns = _infer_introns(result.exons, genome_sequence)

    return result


def _infer_introns(exons, genome_sequence):
    # build intron tuples from the gaps between consecutive exons
    introns = []

    for left, right in zip(exons, exons[1:]):
        start = left[1] + 1
        end = right[0] - 1
        if end >= start:
            introns.append((start, end, genome_sequence[start - 1:end], "."))

    return introns


# VALIDATION

def _preliminary_validate(mergedhit, exonerate_result):
    # cheap check run before the window ladder may stop: catches a prediction far longer
    # than the BLAST hit, or one carrying several homeodomains
    if not exonerate_result.protein_sequence:
        return False

    predicted = exonerate_result.protein_sequence.replace("*", "").replace("-", "")
    original = mergedhit.sequence.replace("*", "").replace("-", "")

    if len(predicted) > len(original) * 3:
        return False

    if len(HELIX3_PATTERN.findall(predicted)) > 2:
        return False

    return True


def _validate_exonerate_result(mergedhit, exonerate_result, args):
    # full check of a prediction against the BLAST hit it came from
    if not exonerate_result.protein_sequence:
        return False

    original = mergedhit.sequence.replace("*", "").replace("-", "")
    predicted = exonerate_result.protein_sequence.replace("*", "").replace("-", "")
    seq_len = len(original)

    # core anchor: 20 residues from the middle of the BLAST hit must still be at least
    # half present in the prediction (the termini are the parts most likely truncated)
    if seq_len > 20:
        core_start = max(0, seq_len // 2 - 10)
        core_end = min(seq_len, seq_len // 2 + 10)
        core = original[core_start:core_end]
    else:
        core = original

    matcher = difflib.SequenceMatcher(None, core, predicted)
    longest = matcher.find_longest_match(0, len(core), 0, len(predicted))

    if longest.size < len(core) * 0.5:
        if args.verbose:
            print(f"      Validation failed: core anchor similarity too low "
                  f"({longest.size}/{len(core)} aa matched)")
        return False

    # length sanity; the 100aa floor stops very short BLAST hits imposing too tight a ceiling
    max_length = max(seq_len * 2, 100)
    if len(predicted) > max_length:
        if args.verbose:
            print(f"      Validation failed: predicted sequence too long "
                  f"({len(predicted)} aa vs max {max_length} aa)")
        return False

    # a single homeodomain per gene, otherwise two genes have been merged
    matches = HELIX3_PATTERN.findall(predicted)
    if len(matches) > 1:
        if args.verbose:
            print(f"      Validation failed: multiple homeodomain signatures "
                  f"({len(matches)} WF[QN]NRR motifs)")
        return False

    # repeated 15-residue fragments mean Exonerate looped over a repeat;
    # only worth checking when the prediction has grown
    if len(predicted) > seq_len * 1.5:
        for i in range(len(predicted) - 15):
            if predicted[i:i + 15] in predicted[i + 15:]:
                if args.verbose:
                    print("      Validation failed: repeated 15-aa fragment detected")
                return False

    return True
