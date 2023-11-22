#!/usr/bin/env python3
"""
download the ENA release from http://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/
"""


import sys
import os
import subprocess
import init
import pathlib
import datetime

COUNT_FILE_PATH = "/home/dbarberis/conta_sequenze.txt"

ena_path = sys.argv[1]

for section in init.ena_sections:
    for type in ["STD", "PAT"]:
        idx = 1
        while True:
            if not pathlib.Path(f"{ena_path}/{type}_{section}_{idx}.fasta").is_file():
                print(f"downloading http://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/{type}_{section}_{idx}.fasta.gz")
                # os.system(f"wget ftp://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/{type}_{section}_{idx}.dat.gz -O {ena_path}/{type}_{section}_{idx}.dat.gz")

                process = subprocess.run(
                    [
                        "wget",
                        f"http://ftp.ebi.ac.uk/pub/databases/ena/sequence/snapshot_latest/std/{type}_{section}_{idx}.fasta.gz",
                        "-O",
                        f"{ena_path}/{type}_{section}_{idx}.fasta.gz",
                    ],
                    capture_output=True,
                )

                if (
                    pathlib.Path(f"{ena_path}/{type}_{section}_{idx}.fasta.gz").is_file()
                    and pathlib.Path(f"{ena_path}/{type}_{section}_{idx}.fasta.gz").stat().st_size == 0
                ):
                    break
                print(f"unzipping {type}_{section}_{idx}.fasta.gz ...")
                os.system(f"gunzip {ena_path}/{type}_{section}_{idx}.fasta.gz")
                print("done")
            else:
                print(f"{type}_{section}_{idx}.fasta.gz already present in {ena_path}")
            idx += 1
            if idx >= 1000:
                break

# count number of sequences
os.system(f"rm -f '{COUNT_FILE_PATH}'")
os.system(f"""for f in {ena_path}/*.fasta; do echo -n "$f " >> '{COUNT_FILE_PATH}'; grep -c '>' $f >> '{COUNT_FILE_PATH}' ; done""")

# write done file
with open(f"/home/dbarberis/DOWNLOAD_DONE_{datetime.datetime.now().isoformat()}", "w") as f_out:
    f_out.write("done")
