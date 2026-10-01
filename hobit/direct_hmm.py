# -*- coding: utf-8 -*-
"""
Protein mode for HobIT: scan protein sequences straight against the homeodomain
HMM, bypassing BLAST, clustering and Exonerate entirely. Used when the input is
already a proteome, so there are only domain boundaries to find.

Here define:
    check_hmmer()
    run_hmmsearch(hmm_path, target_fa, out_domtbl, cpu=1)
    parse_domtbl(path)
    run_protein_pipeline(args, protein_file)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
"""

import shutil
import subprocess
from importlib import resources
from pathlib import Path

from Bio import SeqIO

# domains shorter than this are almost always spurious partial matches
MIN_DOMAIN_LENGTH = 20


# FUNCS FOR RUNNING HMMER

def check_hmmer():
    if not shutil.which("hmmsearch"):
        raise RuntimeError(
            "ERROR: 'hmmsearch' not found in PATH. "
            "Please install HMMER3 (http://hmmer.org/)"
        )


def run_hmmsearch(hmm_path, target_fa, out_domtbl, cpu=1):
    proc = subprocess.run(
        ["hmmsearch",
         "--noali",
         "--cpu", str(cpu),
         "--domtblout", str(out_domtbl),
         str(hmm_path),
         str(target_fa)],
        capture_output=True, text=True
    )

    if proc.returncode != 0:
        raise RuntimeError(f"hmmsearch failed: {proc.stderr}")


def parse_domtbl(path):
    # domtblout columns: 0 target, 3 query, 12 i-Evalue, 13 score, 17 ali_from, 18 ali_to
    hits = []

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.split()
            if len(cols) < 23:
                continue

            hits.append({
                "target": cols[0],
                "query": cols[3],
                "i_evalue": float(cols[12]),
                "score": float(cols[13]),
                "ali_from": int(cols[17]),
                "ali_to": int(cols[18])
            })

    return hits


# MAIN ENTRY POINT FOR PROTEIN MODE

def run_protein_pipeline(args, protein_file):
    check_hmmer()

    input_path = Path(protein_file)

    # a user-supplied profile wins; otherwise fall back to the bundled HD.hmm
    if args.hmm:
        hmm_path = Path(args.hmm)
    else:
        try:
            hmm_path = resources.files("hobit.data").joinpath("HD.hmm")
            if not hmm_path.is_file():
                raise FileNotFoundError(f"Missing HMM file at {hmm_path}")
        except Exception as e:
            print(f"Warning: package-based HMM loading failed ({e}); trying raw path.")
            hmm_path = Path(__file__).parent.joinpath("data", "HD.hmm")

    if not Path(hmm_path).exists():
        print(f"ERROR: HMM profile not found at {hmm_path}")
        print("Please use --hmm to specify the path to HD.hmm")
        return

    base_name = input_path.stem
    out_dir = Path(args.out) / base_name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nProcessing protein file: {base_name}")
    print(f"  - Input: {input_path}")
    print(f"  - HMM: {hmm_path}")
    print(f"  - Output: {out_dir}")

    domtbl_path = out_dir / f"{base_name}.domtblout"
    run_hmmsearch(hmm_path, input_path, domtbl_path, cpu=args.threads)

    # filter on e-value and domain length, then sort best first
    hits = []
    for h in parse_domtbl(domtbl_path):
        domain_length = h["ali_to"] - h["ali_from"] + 1
        if h["i_evalue"] <= args.max_evalue and domain_length >= MIN_DOMAIN_LENGTH:
            hits.append(h)

    hits.sort(key=lambda h: h["i_evalue"])

    tsv_path = out_dir / f"{base_name}_hits.tsv"
    fasta_path = out_dir / f"{base_name}_domains.fasta"

    # indexed rather than loaded: a proteome is large and only the hits are needed
    seq_index = SeqIO.index(str(input_path), "fasta")

    count = 0
    with open(tsv_path, "w") as tsv, open(fasta_path, "w") as fas:
        tsv.write("Protein_ID\tStart\tEnd\tE-value\tScore\tSequence\n")

        for h in hits:
            protein_id = h["target"]
            if protein_id not in seq_index:
                continue

            # write the domain sub-sequence, not the whole protein
            start, end = h["ali_from"], h["ali_to"]
            full_protein = str(seq_index[protein_id].seq)
            domain_seq = full_protein[start - 1:end]

            evalue_str = f"{h['i_evalue']:.2e}"
            score_str = f"{h['score']:.1f}"

            tsv.write(f"{protein_id}\t{start}\t{end}\t{evalue_str}\t{score_str}\t{domain_seq}\n")
            fas.write(f">{protein_id}_{start}-{end} evalue={evalue_str} score={score_str}\n"
                      f"{domain_seq}\n")
            count += 1

    print(f"  - Done! Found {count} domains.")
    print(f"  - Saved to: {tsv_path}")
