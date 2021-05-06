#!/usr/bin/env python3

from Bio import SeqIO
import sys

handle = open(sys.argv[1], "r")

for record in SeqIO.parse(handle, "fasta"):
    description = record.description
    id, ac, year, date = description.split("|")
    # print(f">{record.id}\n{record.seq}")
    print(f">gi|{id}|gisaid|{ac}| {date}")
    print(f"{record.seq}")

handle.close()

