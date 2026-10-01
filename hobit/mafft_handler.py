# -*- coding: utf-8 -*-
"""
MAFFT alignment refinement (LEGACY, enabled with --mafft). HMM refinement in
structures.py is the default and is better in every respect that matters here;
this module is kept because the original HobIT worked this way and some
comparisons still need it.

Here define:
    check_mafft_installed()
    align_sequences_with_mafft(args, mergedhits)
    _align_family_group(args, family_name, group)

@author: Thomas Barton-Owen <t.bartonowen+hobit@gmail.com>, modified by Yurou Ding <yd54@st-andrews.ac.uk>
"""

import subprocess
import tempfile
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def check_mafft_installed():
    # returns True if mafft is on the PATH
    try:
        return subprocess.run(["mafft", "--version"],
                              capture_output=True, text=True).returncode == 0
    except FileNotFoundError:
        return False


def align_sequences_with_mafft(args, mergedhits):
    # align each family's sequences against each other; families with fewer than
    # two members are skipped, since MAFFT has nothing to align them to
    if not check_mafft_installed():
        if args.verbose:
            print("WARNING: MAFFT not found. Skipping MAFFT refinement.")
        return mergedhits

    if args.verbose:
        print("Refining sequences with MAFFT alignment (legacy method)...")

    # group the hits by family
    families = {}
    for mh in mergedhits:
        family = mh.homeodomain_family_guess
        if isinstance(family, (list, tuple)):
            family = family[0]
        families.setdefault(family, []).append(mh)

    refined_count = 0

    for family, group in families.items():
        if len(group) < 2:
            continue

        try:
            _align_family_group(args, family, group)
            refined_count += sum(1 for mh in group if "MAFFT_REFINED" in mh.notes)
        except (subprocess.SubprocessError, ValueError, OSError) as e:
            if args.verbose:
                print(f"Warning: failed to align family {family}: {e}")

    if args.verbose:
        print(f"    Refined {refined_count} sequences using MAFFT")

    return mergedhits


def _align_family_group(args, family_name, group):
    # build one FASTA record per hit in the family
    records = []

    for idx, mh in enumerate(group):
        sequence = mh.improved_sequence or mh.sequence

        # min() of the location, not the tuple itself: "(120, 480)" would put commas
        # and brackets in the FASTA ID, which MAFFT and BioPython choke on
        record_id = f"seq_{idx}_{mh.scaffold}_{min(mh.scaffold_location_guess)}"

        records.append(SeqRecord(Seq(sequence), id=record_id,
                                 description=f"family={family_name}"))

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as temp_in:
        input_path = temp_in.name
        SeqIO.write(records, temp_in, "fasta")

    output_path = input_path + ".aligned.fa"

    try:
        with open(output_path, "w") as out_file:
            subprocess.run(["mafft", "--auto", "--quiet", input_path],
                           stdout=out_file, check=True)

        alignment = AlignIO.read(output_path, "fasta")

        for mh, aligned_record in zip(group, alignment):
            aligned_seq = str(aligned_record.seq).replace("-", "")

            # the alignment may extend a hit but never trim it
            if len(aligned_seq) >= len(mh.sequence):
                mh.improved_sequence = aligned_seq
                if "MAFFT_REFINED" not in mh.notes:
                    mh.notes.append("MAFFT_REFINED")

    finally:
        Path(input_path).unlink(missing_ok=True)
        Path(output_path).unlink(missing_ok=True)

    return group
