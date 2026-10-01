# -*- coding: utf-8 -*-
"""
Compare HobIT predictions against a reference annotation with gffcompare:
run the comparison, write a readable report, mark each hit as Novel or Known,
swap matched predictions for the reference CDS, and optionally drop hits that
fall inside repeats.

Here define:
    check_gffcompare_installed()
    run_gffcompare(args, hobit_gff_path, reference_gff_path)
    _parse_gffcompare_results(args, output_prefix)
    generate_comparison_report(args, results, mergedhits)
    annotate_mergedhits_with_comparison(args, mergedhits, results)
    filter_by_annotation(args, mergedhits, results, filter_repeats=False)
    replace_with_reference_sequence(args, mergedhits, results, ref_gff_path, genome_dict)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
"""

import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from hobit.bioinformatics_tools import BioinformaticsTools

# gffcompare class codes: '=' 'c' 'k' 'j' mean the prediction matches a known
# transcript, 'u' means intergenic (novel), 'r' means it sits in a repeat
MATCHED_CLASS_CODES = ("=", "c", "k", "j")


def check_gffcompare_installed():
    return shutil.which("gffcompare") is not None


# FUNCS FOR RUNNING AND PARSING GFFCOMPARE

def run_gffcompare(args, hobit_gff_path, reference_gff_path):
    # returns the parsed statistics dict, or None if gffcompare is missing or failed;
    # every caller treats None as "skip the comparison steps"
    if not check_gffcompare_installed():
        if args.verbose:
            print("WARNING: gffcompare not found. Install from:")
            print("    https://ccb.jhu.edu/software/stringtie/gffcompare.shtml")
        return None

    if args.verbose:
        print("\nComparing HobIT predictions against reference annotation...")
        print(f"  HobIT GFF: {hobit_gff_path.name}")
        print(f"  Reference: {reference_gff_path.name}")

    comparison_dir = Path(args.out).joinpath("gffcompare_results")
    comparison_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = comparison_dir.joinpath("hobit_vs_ref")

    try:
        result = subprocess.run(
            ["gffcompare",
             "-r", str(reference_gff_path),
             "-o", str(output_prefix),
             str(hobit_gff_path)],
            capture_output=True, text=True, timeout=300
        )
    except subprocess.TimeoutExpired:
        if args.verbose:
            print("WARNING: gffcompare timed out")
        return None

    if result.returncode != 0:
        if args.verbose:
            print(f"WARNING: gffcompare failed: {result.stderr}")
        return None

    if args.verbose:
        print("  Comparison complete")
        print(f"  Results saved to: {comparison_dir}")

    return _parse_gffcompare_results(args, output_prefix)


def _parse_gffcompare_results(args, output_prefix):
    # .stats is only echoed for the user, .tracking gives the prediction-to-reference
    # pairing, .annotated.gtf gives the class code counts
    results = {
        "total_predictions": 0,
        "matched_transcripts": 0,
        "novel_predictions": 0,
        "repeat_overlaps": 0,
        "class_code_counts": defaultdict(int),
        "detailed_matches": []
    }

    stats_file = Path(str(output_prefix) + ".stats")
    if stats_file.exists() and args.verbose:
        print("\n=== gffcompare Summary ===")
        print(stats_file.read_text())

    # .tracking columns: [0] tracking ID, [1] locus, [2] reference ID, [3] class code, [4] query info
    tracking_file = Path(str(output_prefix) + ".tracking")
    if tracking_file.exists():
        with open(tracking_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) < 5:
                    continue

                # an NCBI reference entry reads "gene_id|transcript_id"; the transcript
                # ID (rna-XM_...) is the one CDS features parent to, so keep that half
                ref_raw = cols[2]
                ref_id = "-"
                if ref_raw != "-":
                    parts = ref_raw.split("|")
                    ref_id = parts[1] if len(parts) > 1 else parts[0]

                # query info reads "q1:transcript_id|gene_id|..."; strip the sample
                # prefix and everything after the first pipe to get mRNA_n
                query_id = cols[4]
                if ":" in query_id:
                    query_id = query_id.split(":")[1]
                if "|" in query_id:
                    query_id = query_id.split("|")[0]

                results["detailed_matches"].append({
                    "tracking_id": cols[0],
                    "locus": cols[1],
                    "ref_gene": ref_id if ref_id != "-" else "Novel",
                    "query": query_id
                })

    annotated_file = Path(str(output_prefix) + ".annotated.gtf")
    if annotated_file.exists():
        with open(annotated_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                cols = line.strip().split("\t")
                if len(cols) < 9:
                    continue

                # class_code is written as: class_code "u";
                class_code = None
                for attr in cols[8].split(";"):
                    attr = attr.strip()
                    if attr.startswith("class_code"):
                        class_code = attr.split('"')[1]
                        break

                if not class_code:
                    continue

                results["class_code_counts"][class_code] += 1

                if class_code in MATCHED_CLASS_CODES:
                    results["matched_transcripts"] += 1
                elif class_code == "u":
                    results["novel_predictions"] += 1
                elif class_code == "r":
                    results["repeat_overlaps"] += 1

    results["total_predictions"] = sum(results["class_code_counts"].values())

    return results


# FUNCS FOR REPORTING

def generate_comparison_report(args, results, mergedhits):
    # write comparison_report.txt and echo a short summary to the terminal
    if not results:
        return

    report_path = Path(args.out).joinpath("gffcompare_results", "comparison_report.txt")

    # "or 1" only guards the percentage arithmetic against an empty comparison
    total = results["total_predictions"] or 1

    with open(report_path, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("HobIT vs Reference Annotation Comparison Report\n")
        f.write("=" * 80 + "\n\n")

        f.write("SUMMARY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total HobIT predictions:        {results['total_predictions']}\n")
        f.write(f"Matched to known transcripts:      {results['matched_transcripts']} "
                f"({100 * results['matched_transcripts'] / total:.1f}%)\n")
        f.write(f"Novel/intergenic predictions:      {results['novel_predictions']} "
                f"({100 * results['novel_predictions'] / total:.1f}%)\n")
        f.write(f"Overlapping repetitive regions:    {results['repeat_overlaps']} "
                f"({100 * results['repeat_overlaps'] / total:.1f}%)\n\n")

        f.write("CLASS CODE BREAKDOWN\n")
        f.write("-" * 80 + "\n")
        f.write("Class codes indicate the relationship to reference annotation:\n\n")

        for code in sorted(results["class_code_counts"]):
            count = results["class_code_counts"][code]
            f.write(f"  [{code}] {count:4d} ({100 * count / total:5.1f}%)\n")

        f.write("\n" + "=" * 80 + "\n")

    if args.verbose:
        print(f"\n  Comparison report saved to: {report_path.name}")
        print("\n=== Quick Summary ===")
        print(f"  Matched transcripts: {results['matched_transcripts']}"
              f"/{results['total_predictions']} "
              f"({100 * results['matched_transcripts'] / total:.1f}%)")
        print(f"  Novel predictions:   {results['novel_predictions']} "
              f"({100 * results['novel_predictions'] / total:.1f}%)")
        if results["repeat_overlaps"] > 0:
            print(f"  Repeat overlaps:     {results['repeat_overlaps']} "
                  f"({100 * results['repeat_overlaps'] / total:.1f}%)")


def annotate_mergedhits_with_comparison(args, mergedhits, results):
    # note on each hit whether it matched a reference gene or is novel
    if not results or not results["detailed_matches"]:
        return mergedhits

    match_lookup = {m["query"]: m["ref_gene"] for m in results["detailed_matches"]}

    # the GFF writes hits as gene_n / mRNA_n in list order, so the list position
    # is what links a MergedHit back to its gffcompare row
    for idx, mh in enumerate(mergedhits):
        # gffcompare may report either ID depending on the reference's format
        ref_gene = match_lookup.get(f"gene_{idx + 1}") or match_lookup.get(f"mRNA_{idx + 1}")

        if ref_gene is None:
            continue

        if ref_gene != "Novel":
            mh.notes.append(f"REF_MATCH:{ref_gene}")
        else:
            mh.notes.append("NOVEL_PREDICTION")

    return mergedhits


# FUNCS FOR FILTERING AND RESCUING PREDICTIONS

def filter_by_annotation(args, mergedhits, results, filter_repeats=False):
    # drop predictions that gffcompare put inside a repeat (class code 'r')
    if not results or not filter_repeats:
        return mergedhits

    output_prefix = Path(args.out).joinpath("gffcompare_results", "hobit_vs_ref")
    annotated_file = Path(str(output_prefix) + ".annotated.gtf")

    if not annotated_file.exists():
        return mergedhits

    repeat_indices = set()

    with open(annotated_file) as f:
        for line in f:
            if 'class_code "r"' not in line:
                continue

            match = re.search(r'transcript_id "([^"]+)"', line)
            if not match:
                continue

            # HobIT writes transcripts as mRNA_n, n counting from 1
            transcript_id = match.group(1)
            if transcript_id.startswith("mRNA_"):
                try:
                    transcript_number = int(transcript_id.split("_")[1])
                    repeat_indices.add(transcript_number - 1)
                except ValueError:
                    pass

    if not repeat_indices:
        return mergedhits

    return [mh for i, mh in enumerate(mergedhits) if i not in repeat_indices]


def replace_with_reference_sequence(args, mergedhits, results, ref_gff_path, genome_dict):
    # where gffcompare paired a prediction with a reference transcript, the curated
    # reference CDS replaces the predicted sequence; the original length is noted first
    if not results or not results.get("detailed_matches"):
        return mergedhits

    if args.verbose:
        print("\nAttempting to rescue sequences using reference annotation matches...")

    # list position -> reference transcript ID, same mRNA_n convention as above
    idx_to_ref = {}
    wanted_refs = set()

    for match in results["detailed_matches"]:
        ref_gene = match["ref_gene"]
        if ref_gene in ("-", "Novel"):
            continue

        try:
            query_number = match["query"].split("_")[1]
            idx = int(query_number) - 1
        except (IndexError, ValueError):
            continue

        idx_to_ref[idx] = ref_gene
        wanted_refs.add(ref_gene)

    if not wanted_refs:
        return mergedhits

    if args.verbose:
        print(f"  Parsing reference GFF for {len(wanted_refs)} matched transcripts...")

    # collect every CDS line belonging to a transcript we care about
    ref_cds = {}

    try:
        with open(ref_gff_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue

                # splitting the attribute field into tokens and intersecting with the
                # wanted set finds the transcript without knowing which key holds the ID
                tokens = set(re.split(r'[;="\s]+', parts[8]))
                found = wanted_refs.intersection(tokens)
                if not found:
                    continue

                ref_id = found.pop()
                if ref_id not in ref_cds:
                    ref_cds[ref_id] = {"coords": [], "strand": parts[6],
                                       "scaffold": parts[0]}

                cds_start = int(parts[3])
                cds_end = int(parts[4])
                ref_cds[ref_id]["coords"].append((cds_start, cds_end))

    except OSError as e:
        if args.verbose:
            print(f"  Error reading reference GFF: {e}")
        return mergedhits

    replaced_count = 0

    for idx, ref_id in idx_to_ref.items():
        if idx >= len(mergedhits) or ref_id not in ref_cds:
            continue

        mh = mergedhits[idx]
        data = ref_cds[ref_id]

        # scaffold names may be decorated differently in the two files
        if data["scaffold"] not in mh.scaffold and mh.scaffold not in data["scaffold"]:
            continue

        record = genome_dict.get(mh.scaffold)
        if not record:
            continue

        # GFF coordinates are 1-based and inclusive; splice in genomic order,
        # then reverse complement the whole thing for a minus-strand gene
        coords = sorted(data["coords"], key=lambda c: c[0])

        cds_seq = ""
        for start, end in coords:
            cds_seq += str(record.seq[start - 1:end])

        if data["strand"] == "-":
            cds_seq = BioinformaticsTools.reverse_complement(cds_seq)

        protein_seq = BioinformaticsTools.translate(cds_seq, 1).strip("*")
        if not protein_seq:
            continue

        mh.notes.append(f"ORIGINAL_PRED_LEN:{len(mh.improved_sequence or '')}")
        mh.improved_sequence = protein_seq
        mh.cds_sequence = cds_seq
        mh.notes.append(f"REPLACED_WITH_REF_SEQ:{ref_id}")
        replaced_count += 1

    if args.verbose:
        print(f"  Replaced {replaced_count} predictions with reference GFF sequences.")

    return mergedhits
