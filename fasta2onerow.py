#!/usr/bin/env python3

from Bio import SeqIO
import sys

handle = open(sys.argv[1], "r")

for record in SeqIO.parse(handle, "fasta"):
    # print(f">{record.id}\n{record.seq}")
    print(f">{record.description}\n{record.seq}")

handle.close()

