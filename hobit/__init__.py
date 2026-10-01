# -*- coding: utf-8 -*-
"""
HobIT - Homeobox Identification Tool.
An annotation-independent pipeline for finding homeobox genes in genomes.
Features include TBLASTN search, strand-aware clustering, HMMER domain refinement,
optional Exonerate gene structures prediction, and batch processing.

Usage:
    python -m hobit.hobit_main -g genome.fa -q queries.fa -o output/ --exonerate

@author: Yurou Ding <yd54@st-andrews.ac.uk>
Adapted from Thomas Barton-Owen.
"""

__version__ = "1.0.0"
__license__ = "None"

from hobit import functions
from hobit import structures
from hobit import bioinformatics_tools
from hobit import mafft_handler
from hobit import exonerate_handler
from hobit import batch_processor
from hobit import outputter
from hobit import direct_hmm
from hobit import gff_comparison

__all__ = [
    'functions',
    'structures',
    'gff_comparison',
    'bioinformatics_tools',
    'mafft_handler',
    'exonerate_handler',
    'batch_processor',
    'outputter',
    'direct_hmm'
]