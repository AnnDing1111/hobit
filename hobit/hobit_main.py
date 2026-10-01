# -*- coding: utf-8 -*-
"""
HobIT - Homeobox Identification Tool. Command line entry point.
Genome mode runs the full pipeline per genome: TBLASTN -> cluster -> enrich ->
(optional Exonerate) -> LPORF + HMM refinement -> length recovery -> output.
Protein mode hands off to direct_hmm and skips all of that.

Here define:
    parse_arguments()
    enforce_minimum_length_smart(args, mergedhits, min_len=DEFAULT_MIN_LENGTH)
    process_single_genome(args, genome_path)
    _compare_against_reference(args, mergedhits, hobit_gff_path, genome_dict)
    check_query_format(query_path)
    main()
    _run_batch(args)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
Adapted from the unpublished HobIT prototype by Thomas Barton-Owen.
"""

import argparse
import multiprocessing
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib import resources
from pathlib import Path

from Bio import SeqIO

from hobit import direct_hmm
from hobit.batch_processor import generate_batch_summary
from hobit.bioinformatics_tools import BioinformaticsTools
from hobit.exonerate_handler import run_exonerate_annotation
from hobit.functions import (blast_executer, enrich_mergedhits, file_sys_prep,
                            get_blast_hits, make_mergedhits, progress,
                            serializer, unserializer)
from hobit.mafft_handler import align_sequences_with_mafft
from hobit.outputter import fasta_outputter, gff_outputter, xlsx_outputter

# homeodomains are 60aa; shorter hits get padded back out from their ORF context
DEFAULT_MIN_LENGTH = 60


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="HobIT: Homeobox Identification Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic genome mode
  python -m hobit.hobit_main -g genome.fa -q queries.fa -o output/

  # Protein mode
  python -m hobit.hobit_main -p proteins.fa -o output/

  # With Exonerate
  python -m hobit.hobit_main -g genome.fa -q queries.fa -o output/ --exonerate --gff
        """
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('-g', '--genome', nargs='+',
                      help='Genome file(s) in FASTA format (genome mode)')
    mode.add_argument('-p', '--protein', nargs='+',
                      help='Protein file(s) in FASTA format (protein mode)')

    parser.add_argument('-o', '--out', required=True,
                        help='Output folder path')
    parser.add_argument('-q', '--query',
                        help='Query pool file (default: data/all_HDs.fa, genome mode only)')

    resume = parser.add_argument_group('Resume options')
    resume.add_argument('-x', '--xml',
                        help='Previous TBLASTN XML file')
    resume.add_argument('-pic', '--pickle',
                        help='Previous serialization file')

    analysis = parser.add_argument_group('Analysis options')
    analysis.add_argument('--hmm',
                          help='HMM profile path (default: lib/HD.hmm)')
    analysis.add_argument('-c', '--cluster', type=int, default=200,
                          help='Clustering distance in bp (default: 200)')
    analysis.add_argument('-w', '--window', type=int, default=5000,
                          help='Window size in bp (default: 5000)')
    analysis.add_argument('--mafft', action='store_true',
                          help='Use MAFFT refinement (legacy, HMM preferred)')
    analysis.add_argument('--exonerate', action='store_true',
                          help='Use Exonerate for gene structure prediction')
    analysis.add_argument('--exonerate-model', default='protein2genome',
                          choices=['protein2genome', 'protein2genome:bestfit'],
                          help='Exonerate alignment model (default: protein2genome)')
    analysis.add_argument('--no-mkblastdb', action='store_true',
                          help='Skip automatic BLAST database creation')

    output = parser.add_argument_group('Output options')
    output.add_argument('--no-xlsx', action='store_true',
                        help='Skip Excel output')
    output.add_argument('--gff', action='store_true',
                        help='Output GFF3 format (requires --exonerate)')
    output.add_argument('--max-evalue', type=float, default=0.1,
                        help='E-value threshold (default: 0.1)')
    output.add_argument('--out-introns', action='store_true',
                        help='Include intron sequences (requires --exonerate)')
    output.add_argument('--out-exons', action='store_true',
                        help='Include exon sequences (requires --exonerate)')
    output.add_argument('--fasta-out', nargs='+',
                        choices=['protein', 'cds', 'gene'],
                        help='FASTA output types')

    general = parser.add_argument_group('General options')
    general.add_argument('-v', '--verbose', action='store_true', default=True,
                         help='Detailed progress messages (default: True)')
    general.add_argument('--quiet', action='store_true',
                         help='Suppress progress messages')
    general.add_argument('--threads', type=int, default=1,
                         help='Number of threads (default: 1)')

    comparison = parser.add_argument_group('Annotation comparison options')
    comparison.add_argument('--compare-annotation', metavar='GFF',
                            help='Compare against reference annotation')
    comparison.add_argument('--filter-repeats', action='store_true',
                            help='Filter repetitive regions (requires --compare-annotation)')

    args = parser.parse_args()

    if args.protein:
        args.mode = 'protein'
        args.input_files = args.protein
    else:
        args.mode = 'genome'
        args.input_files = args.genome

        # locate the bundled query pool through the package, falling back to a
        # __file__-relative path when running from an uninstalled source tree
        if not args.query:
            try:
                args.query = str(resources.files('hobit.data').joinpath('all_HDs.fa'))
            except Exception:
                args.query = str(Path(__file__).parent / 'data' / 'all_HDs.fa')

    if args.quiet:
        args.verbose = False

    if (args.gff or args.out_introns or args.out_exons) and not args.exonerate:
        parser.error("--gff, --out-introns, and --out-exons require --exonerate")

    return args


# FUNCS FOR POST-PROCESSING THE HITS

def enforce_minimum_length_smart(args, mergedhits, min_len=DEFAULT_MIN_LENGTH):
    # pad short domains back out to min_len using their surrounding ORF context
    if args.verbose:
        print(f"\nRunning length recovery (target: {min_len}aa)...")

    recovered_count = 0
    skipped_count = 0

    for mh in mergedhits:
        current_seq = mh.improved_sequence or mh.sequence

        if not current_seq or len(current_seq) >= min_len:
            continue

        # the LPORF if there is one, otherwise the raw six-frame translation
        context = mh.lporf or mh.genome_translate

        # usually means process_genome_seq() never found a frame for this hit
        if not context or current_seq not in context:
            mh.notes.append("LENGTH_FIX_FAILED_NO_CONTEXT")
            skipped_count += 1
            continue

        start = context.find(current_seq)
        end = start + len(current_seq)

        # grow the domain symmetrically in both directions
        deficit = min_len - len(current_seq)
        new_start = start - deficit // 2
        new_end = end + (deficit - deficit // 2)

        # push the overflow to the opposite side rather than just clamping,
        # otherwise a domain near either edge comes back short for no reason
        if new_start < 0:
            new_end -= new_start
            new_start = 0

        if new_end > len(context):
            new_start = max(0, new_start - (new_end - len(context)))
            new_end = len(context)

        new_seq = context[new_start:new_end]

        if len(new_seq) > len(current_seq):
            mh.improved_sequence = new_seq
            mh.notes.append(f"SMART_EXTEND: {len(current_seq)}aa->{len(new_seq)}aa")
            recovered_count += 1

            if len(new_seq) < min_len:
                mh.notes.append("ORF_TOO_SHORT_FOR_60AA")

    if args.verbose:
        print(f"    Recovered {recovered_count} short domains")
        if skipped_count > 0:
            print(f"    Skipped {skipped_count} (no valid context)")

    return mergedhits


# THE GENOME PIPELINE

def process_single_genome(args, genome_path):
    genome_name = Path(genome_path).stem

    if args.verbose:
        print(f"\n{'=' * 80}")
        print(f"Processing: {genome_name}")
        print(f"{'=' * 80}\n")

    # a private copy of args per genome, so in batch mode each genome gets its own
    # output folder without the parallel workers treading on each other
    genome_args = argparse.Namespace(**vars(args))
    genome_args.out = str(Path(args.out) / genome_name)
    genome_args.genome = genome_path
    Path(genome_args.out).mkdir(parents=True, exist_ok=True)

    if not args.pickle:
        file_sys_prep(genome_args)

    # keyed on sanitised IDs to match the names BLAST reports; needed even when
    # resuming, for the Exonerate and gffcompare steps
    if args.verbose:
        print(f"Loading genome {genome_path}...")
    try:
        genome_records = SeqIO.parse(genome_path, "fasta")
        genome_dict = SeqIO.to_dict(
            genome_records,
            key_function=lambda rec: BioinformaticsTools.sanitise(rec.id)
        )
    except (OSError, ValueError) as e:
        print(f"FATAL: Error loading genome: {e}")
        sys.exit(1)

    if args.verbose:
        print(f"    Loaded {len(genome_dict)} scaffolds.")

    if args.pickle:
        # resuming: the hits already exist, so jump straight to the outputs
        if args.verbose:
            print("Resuming from serialized data...")
        mergedhits = unserializer(genome_args, args.pickle)

    else:
        progress(genome_args, 5)

        if not args.xml:
            blast_executer(genome_args)
            progress(genome_args, 50)

        hits = get_blast_hits(genome_args)
        progress(genome_args, 60)

        mergedhits = make_mergedhits(genome_args, hits)
        progress(genome_args, 70)

        mergedhits = enrich_mergedhits(genome_args, mergedhits, genome_dict)
        progress(genome_args, 80)

        if genome_args.exonerate:
            # Exonerate is the slow step, so poor hits are filtered out first;
            # the skipped ones rejoin the list afterwards with their BLAST sequences
            if args.verbose:
                print(f"Pre-filtering {len(mergedhits)} hits for Exonerate...")

            hits_to_run = [mh for mh in mergedhits if mh.best_evalue <= genome_args.max_evalue]
            hits_to_skip = [mh for mh in mergedhits if mh.best_evalue > genome_args.max_evalue]

            for mh in hits_to_skip:
                mh.notes.append("EXONERATE_SKIPPED_EVALUE")

            if args.verbose:
                print(f"    Sending {len(hits_to_run)} high-quality hits to Exonerate.")
                print(f"    Skipping {len(hits_to_skip)} low-quality hits.")

            if hits_to_run:
                mergedhits = run_exonerate_annotation(genome_args, hits_to_run, genome_dict) + hits_to_skip
            elif args.verbose:
                print("No hits met the threshold for Exonerate.")

        progress(genome_args, 90)

        # LPORF first, then HMM: improve_homeodomain_sequence() scans every candidate
        # (BLAST hit, LPORF, Exonerate protein), so it must run last, and always runs
        # because it decides where the domain actually starts and ends
        if args.verbose:
            print("Running LPORF generation and HMM refinement...")

        for mh in mergedhits:
            mh.process_genome_seq(genome_args)
            mh.improve_homeodomain_sequence(genome_args)

        mergedhits = enforce_minimum_length_smart(genome_args, mergedhits)

        if args.mafft:
            mergedhits = align_sequences_with_mafft(genome_args, mergedhits)

        # last, once improved_sequence has stopped changing: annotate() flags
        # problems in the sequence that will actually be written out
        for mh in mergedhits:
            mh.annotate()

        serializer(genome_args, mergedhits, "pipeline_state.pkl")

    progress(genome_args, 95)

    if not args.no_xlsx:
        if args.verbose:
            print("Generating Excel output...")
        xlsx_outputter(genome_args, mergedhits)

    if args.gff and args.exonerate:
        if args.verbose:
            print("Generating GFF3 output...")
        hobit_gff_path = gff_outputter(genome_args, mergedhits)

        if genome_args.compare_annotation and hobit_gff_path:
            mergedhits = _compare_against_reference(
                genome_args, mergedhits, hobit_gff_path, genome_dict
            )

    if genome_args.fasta_out:
        if args.verbose:
            print("Generating FASTA outputs...")
        for output_type in genome_args.fasta_out:
            if output_type in ('cds', 'gene') and not genome_args.exonerate:
                if args.verbose:
                    print(f"    Skipping '{output_type}' FASTA (requires --exonerate)")
                continue
            fasta_outputter(genome_args, mergedhits, output_type=output_type)

    progress(genome_args, 100)

    if args.verbose:
        print(f"\nCompleted: {genome_name}\n")

    return mergedhits


def _compare_against_reference(args, mergedhits, hobit_gff_path, genome_dict):
    # gffcompare is optional, so a run without --compare-annotation should not need the module at all
    from hobit import gff_comparison

    ref_path = Path(args.compare_annotation)
    if not ref_path.exists():
        return mergedhits

    if args.verbose:
        print("Running annotation comparison...")

    results = gff_comparison.run_gffcompare(args, hobit_gff_path, ref_path)
    if not results:
        return mergedhits

    gff_comparison.generate_comparison_report(args, results, mergedhits)
    mergedhits = gff_comparison.annotate_mergedhits_with_comparison(args, mergedhits, results)
    mergedhits = gff_comparison.replace_with_reference_sequence(
        args, mergedhits, results, ref_path, genome_dict
    )

    if args.filter_repeats:
        original_count = len(mergedhits)
        mergedhits = gff_comparison.filter_by_annotation(
            args, mergedhits, results, filter_repeats=True
        )
        if args.verbose:
            print(f"  Filtered {original_count - len(mergedhits)} repeat-overlapping predictions")

    return mergedhits


# ENTRY POINT

def check_query_format(query_path):
    # True if any query header carries the "ID|Family|Class" pipes
    try:
        return any('|' in rec.id for rec in SeqIO.parse(query_path, "fasta"))
    except (OSError, ValueError):
        return False


def main():
    args = parse_arguments()

    if args.verbose:
        print("\n" + "=" * 80)
        print("HobIT - Homeobox Identification Tool")
        print("=" * 80 + "\n")

    if args.mode == 'protein':
        for input_file in args.input_files:
            direct_hmm.run_protein_pipeline(args, input_file)
        if args.verbose:
            print("\nAll protein files processed.")
        return

    # without the pipes every hit comes back "Unknown": worth a warning, not a stop
    if not check_query_format(args.query):
        print("=" * 80)
        print(f"WARNING: Query file '{Path(args.query).name}' lacks '|' delimiters")
        print("Expected format: 'ID|Family|Class' (e.g., 'HOX1_HUMAN|Hox1|ANTP')")
        print("Family/Class identification may fail.")
        print("Continuing in 5 seconds...")
        print("=" * 80)
        time.sleep(5)

    if len(args.input_files) == 1:
        process_single_genome(args, args.input_files[0])
    else:
        _run_batch(args)

    if args.verbose:
        print("HobIT analysis complete!\n")


def _run_batch(args):
    # processes rather than threads: each genome loads its own scaffold dict and
    # spends most of its time in BLAST and Exonerate, so there is nothing to share
    if args.verbose:
        print(f"Batch processing {len(args.input_files)} genomes...")
        print(f"Using {args.threads} workers\n")

    # generate_batch_summary() writes into args.out, which otherwise only gets
    # created as a side effect of the per-genome subfolders
    Path(args.out).mkdir(parents=True, exist_ok=True)

    results = {}

    ##
    with ProcessPoolExecutor(max_workers=args.threads,
                             mp_context=multiprocessing.get_context('spawn')) as executor:
        futures = {executor.submit(process_single_genome, args, genome_path): genome_path
                   for genome_path in args.input_files}

        for count, future in enumerate(as_completed(futures), start=1):
            genome_name = Path(futures[future]).stem

            if args.verbose:
                print(f"\n--- Batch progress: {count}/{len(args.input_files)} ---")

            # a genome that raises is recorded as failed and the rest carry on
            try:
                mergedhits = future.result()
                results[genome_name] = {'status': 'success',
                                        'mergedhits': mergedhits,
                                        'count': len(mergedhits)}
                if args.verbose:
                    print(f"Successfully processed: {genome_name}")

            except Exception as e:
                print(f"ERROR processing {genome_name}: {e}")
                traceback.print_exc()
                results[genome_name] = {'status': 'failed',
                                        'error': str(e),
                                        'count': 0}

    generate_batch_summary(args, results)

    return results


if __name__ == "__main__":
    main()
