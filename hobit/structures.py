# -*- coding: utf-8 -*-
"""
Data structures for HobIT. Hit is one BLAST HSP; MergedHit is a cluster of
nearby same-strand Hits, i.e. one candidate gene, and also owns the two
refinement steps that work on a single locus (LPORF and HMM refinement).

Here define:
    check_hmmer_installed()
    missing_hmmpress_files(hmm_path)
    _report_once(key, message)
    get_hmm_db_path()
    class Hit
        __init__(query, query_location, score, evalue, align_list, scaffold,
                 scaffold_location, frame)
    class MergedHit
        __init__(hits)
        annotate()
        cull_hits(num)
        attach_genome_sequence(sequence, offset=0)
        attach_best_query_sequence(sequence)
        process_genome_seq(args)
        improve_homeodomain_sequence(args)
        __list_agree_counter(labels)
        __get_family_from_query(query)
        __get_id_consensus(queries, agree=False)

@author: Yurou Ding <yd54@st-andrews.ac.uk>
Adapted from the unpublished HobIT prototype by Thomas Barton-Owen.
"""

import shutil
import subprocess
import tempfile
from collections import Counter
from importlib import resources
from pathlib import Path

from hobit.bioinformatics_tools import BioinformaticsTools

# hmmscan reads the binary indexes hmmpress writes beside the profile, not the .hmm
# itself (hmmsearch does not, which is why protein mode still works without them)
HMMPRESS_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")

# environment-level problems are the same for every locus, so report them once
_reported_problems = set()


# FUNCS FOR LOCATING HMMER AND ITS PROFILE

def check_hmmer_installed():
    return shutil.which("hmmscan") is not None


def missing_hmmpress_files(hmm_path):
    # the hmmpress suffixes absent beside hmm_path; empty list if it has been pressed
    hmm_path = Path(str(hmm_path))

    missing = []
    for suffix in HMMPRESS_SUFFIXES:
        index_path = hmm_path.with_name(hmm_path.name + suffix)
        if not index_path.exists():
            missing.append(suffix)
    return missing


def _report_once(key, message):
    if key not in _reported_problems:
        _reported_problems.add(key)
        print(message)


def get_hmm_db_path():
    # importlib.resources finds HD.hmm however the package was installed;
    # the __file__-relative fallback covers an uninstalled source tree
    try:
        hmm_path = resources.files("hobit.data").joinpath("HD.hmm")
        if not hmm_path.is_file():
            raise FileNotFoundError(f"HMM file missing at expected path: {hmm_path}")
        return hmm_path

    except Exception as e:
        print(f"Warning: package-based HMM loading failed ({e}); trying raw path.")
        return Path(__file__).parent.joinpath("data", "HD.hmm")


# THE HIT CLASS: ONE BLAST HSP

class Hit:

    # __slots__ keeps the per-hit footprint small; genome-scale runs make
    # hundreds of thousands of these
    __slots__ = ["query", "query_location", "score", "evalue", "align_list",
                 "scaffold", "scaffold_location", "frame"]

    def __init__(self, query, query_location, score, evalue, align_list, scaffold,
                 scaffold_location, frame):
        self.query = query
        self.query_location = query_location
        self.score = score
        self.evalue = evalue
        self.align_list = align_list                # [query_align, match_string, sbjct_align]
        self.scaffold = scaffold
        self.scaffold_location = scaffold_location  # as BLAST gave it: start > end on minus strand
        self.frame = frame                          # -3..+3, never 0; the sign is the strand


# THE MERGEDHIT CLASS: A CLUSTER OF HITS = ONE CANDIDATE GENE

class MergedHit:

    def __init__(self, hits):
        # sorted by e-value, so hits[0] is the representative of the cluster
        self.hits = sorted(hits, key=lambda h: h.evalue)

        # copy the best hit's data onto the MergedHit
        # (scaffold_location_guess may be (high, low) on the minus strand, use min()/max())
        best = self.hits[0]
        self.scaffold = best.scaffold
        self.scaffold_location_guess = best.scaffold_location
        self.frame = best.frame
        self.best_align_list = best.align_list
        self.best_align_query_location = best.query_location
        self.best_evalue = best.evalue
        self.best_score = best.score
        self.best_query = best.query
        self.best_query_sequence = None
        self.query_list = [hit.query for hit in self.hits]

        self.sequence = self.best_align_list[2].replace("-", "")
        self.improved_sequence = self.sequence
        self.notes = []
        self.flag = ""

        # family/class voted across the whole cluster:
        # agree=False gives the winning labels, agree=True gives their vote counts
        self.homeodomain_class_guess, self.homeodomain_family_guess = \
            self.__get_id_consensus(self.query_list, agree=False)
        self.homeodomain_class_guess_agree, self.homeodomain_family_guess_agree = \
            self.__get_id_consensus(self.query_list, agree=True)

        # sequence slots filled in later by the enrichment / annotation steps;
        # genome_sequence is always the forward strand, offset is its 0-based start
        self.genome_sequence = None
        self.genome_window_offset = 0
        self.genome_translate = None
        self.lporf = None

        self.exonerate_result = None
        self.exon_count = 0
        self.intron_count = 0
        self.cds_sequence = ""
        self.full_gene_sequence = ""

        # chimera check: a second high-confidence hit voting for a different family
        # means the cluster probably spans two adjacent genes
        self.is_chimeric = False
        for hit in self.hits[1:]:
            if hit.evalue < 1e-10:
                family = self.__get_family_from_query(hit.query)
                if family != "Unknown" and family != self.homeodomain_family_guess:
                    self.is_chimeric = True
                    self.notes.append(
                        f"WARNING: CHIMERIC_HIT (possible merge: "
                        f"{self.homeodomain_family_guess} + {family})"
                    )
                    break

    def annotate(self):
        # call LAST, after HMM refinement and length recovery: flags stops in the sequence
        # that is actually written out, not in the raw TBLASTN match;
        # a trailing stop is the normal end of a CDS, so only internal ones count
        output_sequence = (self.improved_sequence or self.sequence).rstrip("*")

        # the "not in notes" checks keep this idempotent for resumed runs
        if "*" in output_sequence and "SUBJECT MATCH CONTAINS STOP CODONS" not in self.notes:
            self.notes.append("SUBJECT MATCH CONTAINS STOP CODONS")

        # gaps are a property of the BLAST alignment itself
        if "-" in self.best_align_list[2] and "SUBJECT MATCH CONTAINS GAPS" not in self.notes:
            self.notes.append("SUBJECT MATCH CONTAINS GAPS")

    def cull_hits(self, num):
        # keep only the top num hits in the cluster; num=0 keeps all
        if num > 0:
            self.hits = self.hits[:num]

    def attach_genome_sequence(self, sequence, offset=0):
        # always the raw forward strand: strand is handled downstream
        # by process_genome_seq() (six frames) and by Exonerate itself
        self.genome_sequence = sequence
        self.genome_window_offset = offset

    def attach_best_query_sequence(self, sequence):
        self.best_query_sequence = sequence

    def process_genome_seq(self, args):
        # the anchor is the middle third of the BLAST match, capped at 25 residues
        # (the termini are the parts most likely to be truncated)
        seq_len = len(self.sequence)
        start_anchor = seq_len // 3
        end_anchor = min((seq_len * 2) // 3, start_anchor + 25)

        anchor = self.sequence[start_anchor:end_anchor] or self.sequence
        anchor = anchor.replace("X", "").replace("-", "")
        clean_anchor = anchor.replace("*", "")

        # the frame whose translation contains the anchor is the gene's frame;
        # translate() reverse complements itself for the negative frames
        for frame in (1, 2, 3, -1, -2, -3):
            translation = BioinformaticsTools.translate(self.genome_sequence, frame)

            if clean_anchor and clean_anchor in translation.replace("*", ""):
                self.genome_translate = translation
                break

        if not self.genome_translate:
            self.notes.append("UNEXPECTED PROBLEM FINDING GENOME TRANSLATION")
            self.notes.append("NO_LPORF_FINAL")
            if args.verbose:
                print(f"Unexpected problem finding genome translation for hit at "
                      f"{self.scaffold}:{min(self.scaffold_location_guess)}")
            return

        # the LPORF is the stop-delimited sub-ORF that contains the anchor
        for orf in self.genome_translate.split("*"):
            if clean_anchor and clean_anchor in orf:
                self.lporf = orf
                break

        if not self.lporf:
            # introns or a pseudogene put stops inside the domain: keep the whole
            # translation anyway, a context with stops is better than none
            if clean_anchor in self.genome_translate.replace("*", ""):
                self.lporf = self.genome_translate
                self.notes.append("RESCUED_LPORF_WITH_STOPS")
                if "*" in self.sequence:
                    self.notes.append("BLAST_HIT_CONTAINED_STOPS")
            else:
                self.notes.append("NO_LPORF_FOUND_ANCHOR_MISMATCH")
                self.notes.append("NO_LPORF_FINAL")

    def improve_homeodomain_sequence(self, args):
        # every protein available for this locus (BLAST match, LPORF, Exonerate prediction)
        # is scanned against HD.hmm, and the domain with the best i-Evalue wins
        hmm_db_path = get_hmm_db_path()
        if not hmm_db_path.is_file():
            self.notes.append("IMPROVE_FAILED_HMM_DB_NOT_FOUND")
            _report_once("db_missing",
                         f"    Skipping HMM refinement: {hmm_db_path} not found.")
            return

        if not check_hmmer_installed():
            self.notes.append("IMPROVE_FAILED_HMMER_NOT_FOUND")
            _report_once("hmmer_missing",
                         "    Skipping HMM refinement: hmmscan not found in PATH.")
            return

        ###
        missing = missing_hmmpress_files(hmm_db_path)
        if missing:
            self.notes.append("IMPROVE_FAILED_HMM_NOT_PRESSED")
            _report_once("not_pressed",
                         f"    Skipping HMM refinement: {hmm_db_path} has not been "
                         f"pressed (missing {', '.join(missing)}).\n"
                         f"    Fix with:  hmmpress {hmm_db_path}")
            return

        # args.mafft is only reported, never overwritten: args is shared across every MergedHit
        if args.mafft and "HMM_USED_INSTEAD_OF_MAFFT" not in self.notes:
            if args.verbose:
                print("    (Note: using HMM refinement instead of legacy --mafft)")
            self.notes.append("HMM_USED_INSTEAD_OF_MAFFT")

        # candidates shorter than 20aa are too short for a meaningful domain hit
        candidates = {}
        if self.sequence and len(self.sequence) >= 20:
            candidates["blast"] = self.sequence.replace("*", "").replace("-", "")
        if self.lporf and len(self.lporf) >= 20:
            candidates["lporf"] = self.lporf

        exonerate_protein = self.exonerate_result.protein_sequence if self.exonerate_result else None
        if exonerate_protein and len(exonerate_protein) >= 20:
            candidates["exonerate"] = exonerate_protein

        if not candidates:
            return

        best_seq = None
        best_evalue = 1e5
        best_source = None
        temp_files = []

        try:
            for source, protein_seq in candidates.items():
                with tempfile.NamedTemporaryFile(mode="w", delete=False,
                                                 suffix=".fa") as protein_file:
                    protein_file.write(f">{source}\n{protein_seq}\n")
                    protein_path = protein_file.name

                domtbl_path = protein_path + ".domtbl"
                temp_files.extend([protein_path, domtbl_path])

                subprocess.run(
                    ["hmmscan", "--domtblout", domtbl_path, "--noali",
                     "-E", "1e-5", "--cpu", "1", str(hmm_db_path), protein_path],
                    check=True, capture_output=True, text=True, timeout=120
                )

                if not Path(domtbl_path).exists():
                    continue

                with open(domtbl_path) as f:
                    for line in f:
                        if line.startswith("#"):
                            continue

                        cols = line.split()
                        if len(cols) < 22:
                            continue

                        # domtblout columns (0-based): 12 = i-Evalue,
                        # 17 = ali_from, 18 = ali_to (both 1-based on the target)
                        domain_evalue = float(cols[12])
                        ali_from = int(cols[17])
                        ali_to = int(cols[18])

                        if domain_evalue < best_evalue:
                            best_evalue = domain_evalue
                            best_source = source
                            best_seq = protein_seq[ali_from - 1:ali_to]

            # accept on e-value alone, without comparing lengths: a correctly
            # bounded short domain beats a longer one with ragged edges
            if best_seq and best_evalue < 1e-3:
                self.improved_sequence = best_seq

                # drop any earlier refinement note so only the current one stands
                self.notes = [n for n in self.notes
                              if "IMPROVED_FROM_LPORF" not in n and "HMM_REFINED" not in n]
                self.notes.append(f"HMM_REFINED (source: {best_source}, e-val: {best_evalue:.2e})")

        except subprocess.CalledProcessError as e:
            # hmmscan explains itself on stderr and exits 1; without this
            # the user only gets the exit status, once per locus
            lines = (e.stderr or "").strip().splitlines()
            reason = lines[-1].strip() if lines else f"exit status {e.returncode}"
            self.notes.append("IMPROVE_FAILED_HMMER_ERROR")
            _report_once(f"hmmscan:{reason}",
                         f"    Warning: hmmscan failed ({reason}).\n"
                         f"    Command: {' '.join(e.cmd)}")

        except (subprocess.SubprocessError, OSError, ValueError, IndexError) as e:
            self.notes.append("IMPROVE_FAILED_HMMER_ERROR")
            if args.verbose:
                print(f"    Warning: HMMER refinement failed for {self.scaffold}: {e}")

        finally:
            for path in temp_files:
                Path(path).unlink(missing_ok=True)

    # PRIVATE HELPERS FOR THE FAMILY/CLASS VOTE

    def __list_agree_counter(self, labels):
        # plurality vote: (most common label, its vote count); ("Unknown", 0) if empty
        if not labels:
            return "Unknown", 0

        return Counter(labels).most_common(1)[0]

    def __get_family_from_query(self, query):
        # pull the family out of an "ID|Family|Class" header, degrading gracefully
        parts = str(query).split("|")

        if len(parts) >= 3:
            return parts[-2]      # ID|Family|Class
        elif len(parts) == 2:
            return parts[-1]      # ID|Family
        else:
            return parts[0]       # bare ID, nothing to extract

    def __get_id_consensus(self, queries, agree=False):
        family_list = []
        class_list = []

        for query in queries:
            parts = str(query).split("|")

            if len(parts) >= 3:
                family_list.append(parts[-2])
                class_list.append(parts[-1])
            elif len(parts) == 2:
                family_list.append(parts[-1])
                class_list.append("Unknown")
            else:
                family_list.append(parts[0])
                class_list.append("Unknown")

        best_class = self.__list_agree_counter(class_list)
        best_family = self.__list_agree_counter(family_list)

        # agree=True returns the vote counts instead of the labels
        if agree:
            return best_class[1], best_family[1]

        return best_class[0], best_family[0]
