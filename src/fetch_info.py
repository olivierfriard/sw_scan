"""
Retrieve sequence info from NCBI with the Entrez utility

usage:
python3 fetch_info.py FILE

FILE is a file with one Accession number by line

Must be used with python >= 3.8
"""

import sys
from Bio import Entrez, GenBank
import urllib

Entrez.email = "Davide.Barberis@elitechgroup.com"

with open(sys.argv[1], "r") as f_in:
    for line in f_in:
        try:
            handle = Entrez.efetch(db="nucleotide", id=line.strip(), rettype="gb", retmode="text")
            record = GenBank.read(handle)
            print(f"{line.strip()}\t{record.accession[0]}\t{record.size}\t{record.definition}")
        except urllib.error.HTTPError:
            print(f"{line.strip()} not found in Nucleotide database")
        except Exception:
            print(f"{line.strip()} ERROR")
