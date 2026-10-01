# -*- coding: utf-8 -*-
"""
batch processing functions and summary report generation.
process_multiple_sequences() provides a sequential alternative to the parallel
batch in hobit_main (useful for tight memory or debugging failures).
generate_batch_summary() writes batch_summary.txt at the end of a multi-genome run.

Here define:
    process_multiple_sequences(args, genome_paths, max_workers=None)
    generate_batch_summary(args, results)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
"""

from pathlib import Path


def process_multiple_sequences(args, genome_paths, max_workers=None):
    # returns {genome_name: {status, mergedhits, count}}, same as the parallel path;
    # max_workers is only accepted for signature compatibility, this runs sequentially
    results = {}

    if args.verbose:
        print(f"Processing {len(genome_paths)} genomes...")

    # imported here because hobit_main imports this module (circular otherwise)
    from hobit.hobit_main import process_single_genome

    for genome_path in genome_paths:
        genome_name = Path(genome_path).stem

        # a genome that raises is recorded as failed and the rest carry on
        try:
            if args.verbose:
                print(f"\n  Processing: {genome_name}")

            mergedhits = process_single_genome(args, genome_path)
            results[genome_name] = {
                'status': 'success',
                'mergedhits': mergedhits,
                'count': len(mergedhits)
            }

        except Exception as e:
            results[genome_name] = {
                'status': 'failed',
                'error': str(e),
                'count': 0
            }
            if args.verbose:
                print(f"  ERROR: {e}")

    return results


def generate_batch_summary(args, results):
    # write batch_summary.txt and return its text
    total_genomes = len(results)
    successful = sum(1 for r in results.values() if r['status'] == 'success')
    total_hits = sum(r['count'] for r in results.values() if r['status'] == 'success')

    lines = [
        "=" * 80,
        "HobIT Batch Processing Summary",
        "=" * 80,
        "",
        f"Total genomes processed: {total_genomes}",
        f"  Successful: {successful}",
        f"  Failed: {total_genomes - successful}",
        f"Total merged hits identified: {total_hits}",
        "",
        "Per-genome results:",
        "-" * 80
    ]

    for genome_name, result in sorted(results.items()):
        if result['status'] == 'success':
            lines.append(f"  {genome_name:40s} : {result['count']:4d} hits")
        else:
            lines.append(f"  {genome_name:40s} : FAILED - {result['error']}")

    lines.extend(["", "=" * 80])
    summary_text = "\n".join(lines)

    (Path(args.out) / "batch_summary.txt").write_text(summary_text)

    if args.verbose:
        print(summary_text)

    return summary_text
