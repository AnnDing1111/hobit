# HobIT — Homeobox Identity Tool

**HobIT** is a comprehensive bioinformatics pipeline for the identification, annotation, and classification of homeobox genes in genome sequences. It integrates TBLASTN-based homology search, HMMER domain boundary refinement, and optional Exonerate gene structure prediction into a single, end-to-end automated workflow, with support for both genome-mode and protein-mode analysis and multi-genome batch processing.

---

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [System Requirements](#system-requirements)
4. [Installation](#installation)
5. [Quick Start](#quick-start)
6. [Directory Structure](#directory-structure)
7. [Usage](#usage)
8. [Input File Formats](#input-file-formats)
9. [Output Files](#output-files)
10. [Parameter Reference](#parameter-reference)
11. [Pipeline Architecture](#pipeline-architecture)
12. [Troubleshooting](#troubleshooting)
13. [Citation](#citation)

---

## Overview

Homeobox genes encode homeodomain transcription factors that play fundamental roles in developmental gene regulation across the Metazoa. Identifying and classifying homeodomain-containing sequences in newly sequenced genomes is a routine yet technically demanding task, requiring integration of sequence homology search, reading-frame–aware translation, and profile-based domain boundary detection.

HobIT automates this process. Starting from a genome assembly (FASTA) and a curated protein query library, HobIT performs TBLASTN search, clusters overlapping hits into candidate loci, extracts genomic context windows, translates all six reading frames, identifies the Longest Possible ORF (LPORF) containing the homeodomain anchor, and refines domain boundaries using HMMER against the Pfam PF00046 homeobox HMM profile. Optionally, Exonerate can be invoked to generate spliced gene structure models, producing GFF3 output suitable for downstream genome annotation pipelines.

---

## Features

**Core pipeline (Genome Mode)**
* TBLASTN-based homology search against nucleotide genome sequences using a curated homeodomain protein query library.
* Automatic BLAST+ database creation; compatible with both BLAST+ v4 and v5 index formats (`.nhr`, `.nsi`, `.nsq`).
* Strand-aware hit clustering: merges overlapping and adjacent HSPs from the same strand into candidate gene loci, preventing incorrect merging of hits from opposite strands.
* Chimeric hit detection: flags loci where high-confidence hits map to multiple homeodomain families.
* Six-frame genomic window translation and LPORF extraction, with context-guided anchor-matching.
* HMMER-based domain boundary refinement (PF00046 profile): scans BLAST match sequence, LPORF, and Exonerate protein independently.
* Smart sequence length recovery: symmetrically expands short domain sequences (< 60 aa) within LPORF or translated genome context.
* Family and class prediction from structured query headers (`ID|Family|Class` format).
* Optional Exonerate gene structure prediction (`protein2genome` model), with E-value pre-filtering.
* Optional annotation comparison via gffcompare, with overlap classification (Novel / Known) and repeat filtering.

**Protein Mode**
* Direct hmmsearch of protein FASTA input against the PF00046 HMM profile, bypassing the TBLASTN pipeline.
* Domain table parsing with E-value and minimum-length filtering.
* Per-protein TSV and FASTA output.

**Batch Processing**
* Multi-genome processing with `ProcessPoolExecutor`-based parallelism.
* Per-genome output subdirectories.
* Pipeline state serialization (pickle) for resume capability.

**Output & Web Interface**
* User-friendly Web GUI for point-and-click analysis.
* Excel workbook (`.xlsx`) with sortable, filterable columns, E-value formatting, and optional CDS/gene sequences.
* GFF3 gene structure annotations with hierarchical feature records (gene → mRNA → exon/CDS/intron).
* FASTA export for protein (refined domain), CDS, or full gene sequences.
* Optional AI-based similarity scoring via ESM2 (`esm_plugin.py`, beta).

---

## System Requirements

**Operating system:** Linux, macOS, or Windows (via WSL2 recommended for bioinformatics tools).  
**Environment Manager:** Conda or Mamba (Miniconda/Anaconda).

HobIT uses a Conda environment to automatically manage all Python dependencies and external bioinformatics tools (BLAST+, HMMER3, Exonerate, MAFFT, gffcompare). **You do not need to install these tools manually.**

---

## Installation

HobIT can be installed and deployed in minutes using Conda and pip.

### 1. Clone the repository

```bash
git clone https://github.com/your-repository/hobit.git
cd hobit/
```

### 2. Create the environment and install

Create the isolated environment containing all bioinformatics dependencies, then install HobIT as a system-wide command-line tool:

```bash
conda env create -f environment.yml
conda activate hobit
pip install .
```

> **Note:** If you plan to use the ESM2 AI scoring plugin (beta), you must install the additional deep learning libraries manually:
> ```bash
> pip install torch transformers scipy pandas numpy
> ```

---

## Quick Start

A minimal "toy dataset" (`tests/toy_genome.fa`, containing a genomic scaffold from the sea anemone *Nematostella vectensis* GCF_932526225.1_jaNemVect1.1_genomic.fa) is provided to verify your installation.

Run the following command from the root directory of the repository to test the pipeline:

```bash
hobit -g tests/toy_genome.fa -o test_output/ --exonerate --fasta-out protein --max-evalue 1e-5
```

This test should complete in about 5 minutes on a standard laptop. A successful run will generate the directory `test_output/` containing `hobit_output.xlsx`and the extracted homeodomain protein sequences in `hobit_protein.fasta`.

---

## Directory Structure

```
HobIT/
├── environment.yml           # Conda dependencies configuration
├── pyproject.toml            # Python packaging configuration
├── tests/                    # Toy datasets for sanity checks
├── hobit/                    # Main package source code
│   ├── hobit_main.py         # CLI entry point
│   ├── web_gui.py            # GUI entry point
│   ├── templates/            # HTML templates for the Web GUI
│   ├── data/                 # Built-in data (automatically loaded)
│   │   ├── all_HDs.fa        # Default homeodomain query library
│   │   ├── HD.hmm            # Pfam homeobox HMM profile
│   │   └── HD.hmm.h3*        # HMMER index files
│   └── ...                   # Core modules
└── README.md
```

---

## Usage

Once installed, HobIT provides two global terminal commands: `hobit` (for CLI processing) and `hobit-gui` (for the web interface).

### Graphical User Interface (Recommended)

Simply type the following command anywhere in your terminal:

```bash
hobit-gui
```

The application will automatically find an available port and open the web interface in your default browser.

### Command-Line Interface (CLI)

Basic genome mode:

```bash
# The default query library and HMM models are loaded automatically!
hobit -g genome.fa -o output/
```

Genome mode with Exonerate and GFF3 output:

```bash
hobit -g genome.fa -o output/ --exonerate --gff --fasta-out protein cds
```

Batch processing (multiple genomes):

```bash
hobit -g genome1.fa genome2.fa genome3.fa -o batch_output/ --threads 4
```

Protein mode (direct HMM scan):

```bash
hobit -p proteins.fa -o output/
```

Resume from BLAST XML or Pickle:

```bash
hobit -g genome.fa -o output/ -x output/genome/blast_output.xml
# OR
hobit -g genome.fa -o output/ --pickle output/genome/pipeline_state.pkl
```

---

## Input File Formats

**Genome file (`-g`):** Standard FASTA format. Scaffold identifiers are sanitized automatically to remove common NCBI prefixes (`ref|`, `lcl|`, `gb|`, `emb|`). Multiple files may be provided for batch processing.

**Query protein library (`-q`):** FASTA format with structured headers. The `|`-delimited header format enables automatic family and class annotation of predicted hits:

```
>UniqueID|FamilyName|ClassName
MTEYKLVVVGAGGVGKSALTIQLIQNHFV...
```

**Protein input (`-p`, protein mode):** Standard FASTA format. Sequences are scanned directly with hmmsearch against PF00046.

**Reference annotation (`--compare-annotation`):** GFF3 or GTF format, compatible with gffcompare.

---

## Output Files

All output files are written to the directory specified with `-o`. In batch mode, each genome receives its own subdirectory.

| File | Description |
|------|-------------|
| `hobit_output.xlsx` | Excel workbook with all merged hits, sorted by E-value |
| `hobit_output.gff3` | GFF3 gene structure annotations (requires `--exonerate --gff`) |
| `hobit_protein.fasta` | Refined homeodomain protein sequences (requires `--fasta-out protein`) |
| `hobit_cds.fasta` | CDS nucleotide sequences (requires `--exonerate --fasta-out cds`) |
| `hobit_gene.fasta` | Full gene sequences including introns (requires `--exonerate --fasta-out gene`) |
| `blast_output.xml` | Raw TBLASTN output (retained for resume capability) |
| `pipeline_state.pkl` | Serialized pipeline state for resume (`--pickle`) |
| `gffcompare_results/` | gffcompare outputs (requires `--compare-annotation`) |

---

## Parameter Reference

### Mode selection (one required)

| Argument | Description |
|----------|-------------|
| `-g` / `--genome FILE [...]` | Genome FASTA file(s) — genome mode |
| `-p` / `--protein FILE [...]` | Protein FASTA file(s) — protein mode |

### Required & Resume

| Argument | Default | Description |
|----------|---------|-------------|
| `-o` / `--out DIR` | — | Output directory (created if absent) |
| `-x` / `--xml FILE` | — | Resume from existing TBLASTN XML output |
| `--pickle FILE` | — | Resume from serialized pipeline state |

### Analysis options

| Argument | Default | Description |
|----------|---------|-------------|
| `-c` / `--cluster INT` | 200 | Clustering distance in bp |
| `-w` / `--window INT` | 5000 | Genomic context window size (±bp) |
| `--exonerate` | off | Enable Exonerate gene structure prediction |
| `--exonerate-model` | `protein2genome` | Exonerate alignment model |
| `--no-mkblastdb` | off | Skip automatic BLAST database creation |
| `--threads INT` | 1 | CPU threads for BLAST and batch processing |


### Output options

| Argument | Default | Description |
|----------|---------|-------------|
| `--gff` | off | Generate GFF3 output (requires `--exonerate`) |
| `--max-evalue FLOAT` | 0.1 | E-value threshold for output filtering |
| `--out-introns` | off | Include full gene sequences |
| `--out-exons` | off | Include CDS sequences |
| `--fasta-out {protein,cds,gene}` | — | FASTA output type(s) |

---

## Pipeline Architecture

The HobIT genome-mode pipeline proceeds through the following stages:

1. Input validation and output directory preparation
2. TBLASTN search (or load from `--xml`)
3. BLAST XML parsing → Hit objects
4. Strand-aware hit clustering → MergedHit objects
5. Genomic window retrieval (±window bp, forward strand)
6. Query sequence retrieval and attachment
7. [Optional] Exonerate gene structure prediction
8. LPORF generation (6-frame translation, anchor matching)
9. HMMER domain boundary refinement (PF00046)
10. Smart sequence length recovery (< 60 aa expansion)
11. Pipeline state serialization
12. Output generation (Excel, GFF3, FASTA)
13. [Optional] Annotation comparison (gffcompare)

---

## Troubleshooting

**`WARNING: Query file lacks '|' delimiters`**  
The query FASTA headers are not in `ID|Family|Class` format. Family and class predictions will be unavailable. Analysis continues after a 5-second pause.

**`UNEXPECTED PROBLEM FINDING GENOME TRANSLATION` in the Notes column**  
The anchor sequence from the BLAST hit could not be located in any of the six reading-frame translations of the genomic window. This typically occurs for heavily intron-interrupted genes, pseudogenes, or loci at scaffold boundaries. The hit is retained with the raw BLAST sequence.

**Exonerate produces no predictions for a locus**  
The locus is retained in output with the BLAST-derived or HMMER-refined sequence. GFF3 output only includes loci with successful Exonerate predictions. Reduce `--max-evalue` to limit Exonerate to higher-confidence hits if runtime is a concern.

---

## Citation

If you use HobIT in your research, please cite:

> *Manuscript in preparation.*

HobIT is provided without warranty. For questions or bug reports, please contact the authors.