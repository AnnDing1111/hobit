# -*- coding: utf-8 -*-
"""
Output generation for HobIT: the Excel table, GFF3 gene structures and FASTA
sequences.

Here define:
    xlsx_outputter(args, mergedhits)
    _column_title_constructor(args)
    _row_data_constructor(args, mergedhit, row)
    _first_of(guess)
    _frame_formatter(frame)
    _gff_line(scaffold, feature, start, end, score, strand, phase, attributes)
    gff_outputter(args, mergedhits)
    fasta_outputter(args, mergedhits, output_type="protein")

@author: Yurou Ding <yd54@st-andrews.ac.uk>
Adapted from the unpublished HobIT prototype by Thomas Barton-Owen.
"""

from pathlib import Path

import xlsxwriter


# FUNCS FOR EXCEL OUTPUT

def xlsx_outputter(args, mergedhits):
    # sort the list of mergedhits by evalue, best first
    mergedhits.sort(key=lambda mh: mh.best_evalue)

    output_path = Path(args.out).joinpath("hobit_output.xlsx")

    # instantiate the workbook & worksheet
    workbook = xlsxwriter.Workbook(output_path)
    worksheet = workbook.add_worksheet("HobIT output")

    header_fmt = workbook.add_format({
        "bold": True, "font_name": "Calibri", "align": "left", "bg_color": "#D3D3D3"
    })
    monospace_fmt = workbook.add_format({"font_name": "Courier New"})
    center_fmt = workbook.add_format({"align": "center"})
    number_fmt = workbook.add_format({"num_format": "0.00E+00"})

    # which columns exist depends on the flags this run was given
    titles = _column_title_constructor(args)
    worksheet.write_row(0, 0, titles)

    row = 1
    for mh in mergedhits:
        row_data = _row_data_constructor(args, mh, row)
        if row_data:                      # None means the hit was filtered out
            worksheet.write_row(row, 0, row_data)
            row += 1

    # column widths, in the same order the titles were built
    col = 0
    for width, fmt in [(12, None),        # Scaffold
                       (10, None),        # Location
                       (4.5, None),       # Frame
                       (12, None),        # Family
                       (10, None),        # Class
                       (10, number_fmt),  # e-value
                       (7, center_fmt),   # Score
                       (5.2, center_fmt)]:  # Length
        worksheet.set_column(col, col, width, fmt)
        col += 1

    if args.exonerate:
        worksheet.set_column(col, col, 6, center_fmt)      # Exons
        col += 1
        worksheet.set_column(col, col, 6, center_fmt)      # Introns
        col += 1

    worksheet.set_column(col, col, 70, monospace_fmt)      # Sequence
    col += 1

    if args.exonerate and args.out_exons:
        worksheet.set_column(col, col, 70, monospace_fmt)  # CDS
        col += 1

    worksheet.set_column(col, col, 70, monospace_fmt)      # ORF
    col += 1

    if args.exonerate and args.out_introns:
        worksheet.set_column(col, col, 70, monospace_fmt)  # Full gene
        col += 1

    worksheet.set_column(col, col, 70, monospace_fmt)      # Genome window
    col += 1
    worksheet.set_column(col, col, 50)                     # Notes

    worksheet.set_row(0, 20, header_fmt)
    worksheet.freeze_panes(1, 0)
    worksheet.autofilter(0, 0, row - 1, len(titles) - 1)

    workbook.close()

    if args.verbose:
        print(f"Exported data to {output_path.name}")
        print(f"    {row - 1} merged hits written")


def _column_title_constructor(args):
    titles = ["Scaffold", "Location", "Fr.", "Family (?)", "Class (?)",
              "e-value", "Score", "Len"]

    if args.exonerate:
        titles.extend(["Exons", "Introns"])

    titles.append("Sequence")

    if args.exonerate and args.out_exons:
        titles.append("CDS Sequence")

    titles.append("ORF")

    if args.exonerate and args.out_introns:
        titles.append("Full Gene Sequence")

    titles.extend(["Nucleotide sequence window", "Notes"])

    return titles


def _row_data_constructor(args, mergedhit, row):
    # e-value is the only filter: hits without an Exonerate prediction are kept on
    # purpose, dropping them cost real genes (Ventl, Cdx) whose domains were fine
    if mergedhit.best_evalue > args.max_evalue:
        return None

    family = _first_of(mergedhit.homeodomain_family_guess)
    class_ = _first_of(mergedhit.homeodomain_class_guess)

    # the Length cell is a live =LEN() formula over the sequence cell, so its column
    # letter depends on the optional exonerate columns sitting in front of it
    seq_col_index = 10 if args.exonerate else 8
    seq_col = chr(ord("A") + seq_col_index)
    length = f"=LEN({seq_col}{row + 1})"

    # the leftmost coordinate is reported on both strands
    data = [mergedhit.scaffold,
            min(mergedhit.scaffold_location_guess),
            _frame_formatter(mergedhit.frame),
            family,
            class_,
            mergedhit.best_evalue,
            mergedhit.best_score,
            length]

    if args.exonerate:
        data.extend([mergedhit.exon_count, mergedhit.intron_count])

    data.append(mergedhit.improved_sequence or mergedhit.sequence)

    if args.exonerate and args.out_exons:
        data.append(mergedhit.cds_sequence)

    data.append(mergedhit.lporf or mergedhit.genome_translate or "")

    if args.exonerate and args.out_introns:
        data.append(mergedhit.full_gene_sequence)

    data.append(mergedhit.genome_sequence or "")
    data.append("; ".join(mergedhit.notes))
    return data


def _first_of(guess):
    # family/class guesses are usually a string but may arrive as a tuple
    if isinstance(guess, (list, tuple)):
        return guess[0]

    return guess


def _frame_formatter(frame):
    # render a reading frame as "+1" / "-3", or "?" if it is out of range
    if 1 <= frame <= 3:
        return f"+{frame}"
    elif -3 <= frame <= -1:
        return str(frame)

    return "?"


# FUNCS FOR GFF3 OUTPUT

def _gff_line(scaffold, feature, start, end, score, strand, phase, attributes):
    # one GFF3 line: seqid, source, type, start, end, score, strand, phase, attributes
    columns = [scaffold, "HobIT", feature, start, end, score, strand, phase, attributes]

    return "\t".join(str(column) for column in columns) + "\n"


def gff_outputter(args, mergedhits):
    # returns the output path, which the caller passes on to gffcompare
    if not args.exonerate:
        if args.verbose:
            print("WARNING: GFF output requires --exonerate flag")
        return None

    output_path = Path(args.out).joinpath("hobit_output.gff3")

    genome = args.genome[0] if isinstance(args.genome, list) else args.genome

    with open(output_path, "w") as gff:
        gff.write("##gff-version 3\n")
        gff.write("# HobIT gene structure predictions\n")
        gff.write(f"# Query: {Path(args.query).name}\n")
        gff.write(f"# Genome: {Path(genome).name}\n")
        gff.write(f"# Max e-value: {args.max_evalue}\n")
        gff.write("#\n")

        gene_id = 1

        for mh in mergedhits:
            if not mh.exonerate_result or not mh.exonerate_result.exons:
                continue
            if mh.best_evalue > args.max_evalue:
                continue

            result = mh.exonerate_result

            # Exonerate coordinates are 1-based and relative to the genome window;
            # adding the window's 0-based start gives the 1-based genomic position
            offset = mh.genome_window_offset

            family = _first_of(mh.homeodomain_family_guess)
            class_ = _first_of(mh.homeodomain_class_guess)

            gene_start = result.gene_start + offset
            gene_end = result.gene_end + offset

            # gene line
            gene_attributes = (f"ID=gene_{gene_id};Name={family};Family={family};"
                               f"Class={class_};evalue={mh.best_evalue:.2e}")
            gff.write(_gff_line(mh.scaffold, "gene", gene_start, gene_end,
                                result.score, result.strand, ".", gene_attributes))

            # mRNA line
            mrna_attributes = f"ID=mRNA_{gene_id};Parent=gene_{gene_id}"
            gff.write(_gff_line(mh.scaffold, "mRNA", gene_start, gene_end,
                                result.score, result.strand, ".", mrna_attributes))

            # exon and CDS come from the same coordinates: protein2genome
            # only ever reports coding exons
            for feature in ("exon", "CDS"):
                prefix = feature.lower()
                for idx, (start, end, _, phase) in enumerate(result.exons, 1):
                    attributes = f"ID={prefix}_{gene_id}_{idx};Parent=mRNA_{gene_id}"
                    gff.write(_gff_line(mh.scaffold, feature, start + offset, end + offset,
                                        ".", result.strand, phase, attributes))

            # intron lines
            for idx, (start, end, _, _) in enumerate(result.introns, 1):
                attributes = f"ID=intron_{gene_id}_{idx};Parent=mRNA_{gene_id}"
                gff.write(_gff_line(mh.scaffold, "intron", start + offset, end + offset,
                                    ".", result.strand, ".", attributes))

            gene_id += 1

    if args.verbose:
        print(f"Exported GFF3 annotations to {output_path.name}")
        print(f"    {gene_id - 1} genes with exon/intron structures")

    return output_path


# FUNCS FOR FASTA OUTPUT

def fasta_outputter(args, mergedhits, output_type="protein"):
    # output_type is 'protein' (the refined domain), 'cds' or 'gene';
    # the latter two only exist when Exonerate ran
    output_path = Path(args.out).joinpath(f"hobit_{output_type}.fasta")

    with open(output_path, "w") as fasta:
        for mh in mergedhits:
            if mh.best_evalue > args.max_evalue:
                continue

            if output_type == "protein":
                sequence = mh.improved_sequence or mh.sequence
            elif output_type == "cds":
                sequence = mh.cds_sequence
            elif output_type == "gene":
                sequence = mh.full_gene_sequence
            else:
                sequence = ""

            # hits with an empty sequence are skipped rather than written as empty records
            if not sequence:
                continue

            family = _first_of(mh.homeodomain_family_guess)
            location = min(mh.scaffold_location_guess)
            frame = _frame_formatter(mh.frame)

            fasta.write(f">{mh.scaffold}_{location} family={family} "
                        f"evalue={mh.best_evalue:.2e} frame={frame}\n")

            # 60 residues per line
            for i in range(0, len(sequence), 60):
                fasta.write(sequence[i:i + 60] + "\n")

    if args.verbose:
        print(f"Exported {output_type} sequences to {output_path.name}")
