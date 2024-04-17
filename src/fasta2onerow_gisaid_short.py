#!/usr/bin/env python3

from Bio import SeqIO
import sys

handle = open(sys.argv[1], "r")

for record in SeqIO.parse(handle, "fasta"):
    description = record.description   # >hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020|2020-04-17

    id, ac, year, date = description.split("|")

    short_ac = ac.replace('EPI_ISL_', '')

    id_split = id.split("/")

    if len(id_split) == 4:
        _, country, code, year = id.split("/")
    elif len(id_split) == 5:
        country = id_split[1] + '_' + id_split[2]
        code = id_split[3] 
        year = id_split[4] 
    else:
        print(id)
        raise

    country_nospace = country.replace(' ', '_')

    # print(f">{record.id}\n{record.seq}")

    print(f">{country_nospace}/{code}|{short_ac}")
    print(f"{record.seq}")

handle.close()

