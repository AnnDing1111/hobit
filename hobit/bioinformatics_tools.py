# -*- coding: utf-8 -*-
"""
Sequence manipulation utilities for HobIT: reverse complement, translation
and ID sanitisation, all as classmethods on one namespace class.

Here define:
    class BioinformaticsTools:
        reverse_complement(sequence)
        sanitise(name)
        translate(sequence, frame=1)
@author: Thomas Barton-Owen <t.bartonowen+hobit@gmail.com>, modified by Yurou Ding <yd54@st-andrews.ac.uk>
"""


class BioinformaticsTools:

    # complement table: DNA, RNA (U -> A) and the IUPAC ambiguity codes
    __complement_table = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "U": "A", "M": "K", "R": "Y", "W": "W",
        "S": "S", "Y": "R", "K": "M", "V": "B",
        "H": "D", "D": "H", "B": "V", "N": "N"
    }

    # standard genetic code, plus the fourfold-degenerate NNN codons whose
    # amino acid is unambiguous
    __codon_table = {
        "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
        "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
        "TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
        "TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",
        "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
        "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
        "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
        "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
        "ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
        "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
        "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
        "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
        "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
        "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
        "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
        "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
        "CTN": "L", "GTN": "V", "ACN": "T", "GCN": "A",
        "CCN": "P", "TCN": "S", "GGN": "G", "CGN": "R"
    }

    @classmethod
    def reverse_complement(cls, sequence):
        # anything not in the complement table becomes N
        sequence = str(sequence).upper()

        return "".join(cls.__complement_table.get(base, "N") for base in reversed(sequence))

    @classmethod
    def sanitise(cls, name):
        # BLAST XML gives the whole description line ("KXJ24917.1 hypothetical protein")
        # while SeqIO gives only "KXJ24917.1"; keeping the first token is what lets the
        # two match, otherwise the query sequence is silently empty and Exonerate is skipped
        name = str(name)

        for prefix in ["ref|", "lcl|", "gb|", "emb|", " No definition line"]:
            name = name.replace(prefix, "")

        tokens = name.rstrip("|").split()

        return tokens[0] if tokens else ""

    @classmethod
    def translate(cls, sequence, frame=1):
        # frame is 1-3 or -1--3; negative frames are reverse complemented first
        assert frame in (-3, -2, -1, 1, 2, 3), "Invalid frame (must be +-1, +-2 or +-3)"

        sequence = str(sequence).upper().replace("U", "T")

        if frame < 0:
            sequence = cls.reverse_complement(sequence)

        offset = abs(frame) - 1  # frame 1 -> offset 0
        sequence = sequence[offset:]

        protein = []
        for i in range(0, len(sequence) - 2, 3):
            # push anything exotic to N, so the lookup can still hit a degenerate NNN entry
            codon = ""
            for base in sequence[i:i + 3]:
                codon += base if base in "ATCGN" else "N"

            # unresolvable codons become X instead of raising, which keeps non-model-species genomes from crashing the pipeline
            protein.append(cls.__codon_table.get(codon, "X"))

        return "".join(protein)
