#!/usr/bin/env python3

"""
convert multi-row fasta fiel in a one-row fasta file

Usage:

python3 fasta2onerow.py MULTI_ROWS_FASTA_FILE > ONE_ROW_FASTA_FILE

"""

from Bio import SeqIO
import sys

handle = open(sys.argv[1], "r")

for record in SeqIO.parse(handle, "fasta"):
    # print(f">{record.id}\n{record.seq}")
    print(f">{record.description}\n{record.seq}")

handle.close()

