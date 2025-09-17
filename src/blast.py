#!/opt/python-3.8.11/bin/python3
"""
Wrapper for BLAST

The BLAST TSV results file is converted to a sqlite database
"""

import sys
import subprocess
import argparse
import pathlib as pl
import os
import uuid

WS_DEFAULT = 10
MAX_TARGET_SEQS = 10_000_000
E_VALUE = 10_000
BLAST_COMMAND = "blastn"
TMP_DIR = "/data/tmp"

"""
available fields are:
id  description frame identity score align_length  target_length aligned_query_sequence aligned_target_sequence  query_begin  query_end  target_begin  target_end_optimal
"""

PARAMETERS = "sseqid stitle sframe pident score length slen qseq sseq qstart qend sstart send evalue bitscore"


__version__ = "4"
__version_date__ = "2022-05-19"

parser = argparse.ArgumentParser(
    prog="BLAST wrapper",
    description="BLAST",
    usage=("blast.py -q QUERY_PATH -t TARGET_PATH -c N_CORES -o OUTPUT_PATH"),
)
parser.add_argument(
    "-q",
    "--query",
    action="store",
    dest="query",
    type=str,
    help="Path of the query file (FASTA format)",
)
parser.add_argument(
    "-t", "--target", action="store", dest="target", type=str, help="division"
)

parser.add_argument(
    "-w",
    "--ws",
    action="store",
    dest="ws",
    default=WS_DEFAULT,
    type=int,
    help="word size",
)

parser.add_argument(
    "-c",
    "--cpu",
    action="store",
    dest="cpu",
    default=8,
    type=int,
    help="Set number of CPU/cores to use (default 8)",
)

parser.add_argument(
    "-o",
    "--output",
    action="store",
    dest="output",
    type=str,
    help="Set path for the output file",
)


parser.add_argument(
    "--temp",
    action="store",
    dest="temp",
    default=TMP_DIR,
    type=str,
    help="Set path for the temporary directory",
)

parser.add_argument(
    "-v",
    "--version",
    action="version",
    version=f"%(prog)s v.{__version__} {__version_date__} (c) Olivier Friard 2021",
    help="Display the help",
)

args = parser.parse_args()

if not args.query:
    print("The query is missing.")
    sys.exit()


if not args.target:
    print("The database is missing.")
    sys.exit()


if not args.output:
    print("The output file is missing.")
    sys.exit()

if args.temp != TMP_DIR:
    TMP_DIR = args.temp

# check if tmp dir exists
if not pl.Path(TMP_DIR).is_dir():
    print(f"The temporary directory {TMP_DIR} does not exist")
    sys.exit()

if args.ws:
    word_size = args.ws
else:
    word_size = WS_DEFAULT

tmp_output_name = str(pl.Path(TMP_DIR) / str(uuid.uuid4()))

command = [
    BLAST_COMMAND,
    "-db",
    args.target,
    "-query",
    args.query,
    "-evalue",
    str(E_VALUE),
    "-max_target_seqs",
    str(MAX_TARGET_SEQS),
    "-word_size",
    str(word_size),
    "-out",
    tmp_output_name,
    "-num_threads",
    str(args.cpu),
    "-outfmt",
    f"6 {PARAMETERS}",
]

print(command)

subprocess.run(command)

# with open(TMP_OUTPUT_NAME + ".nl", "w") as f_out:
#    subprocess.run(["nl", TMP_OUTPUT_NAME], stdout=f_out)
# os.remove(TMP_OUTPUT_NAME)


create_table_cmd = (
    "CREATE TABLE sequences ("
    # "id INTEGER PRIMARY KEY AUTOINCREMENT UNIQUE NOT NULL, "
    "id text,"
    "description text,"
    "frame int,"
    "identity float,"
    "score float,"
    "align_length int,"
    "target_length int,"
    "aligned_query_sequence text,"
    "aligned_target_sequence text,"
    "query_begin int,"
    "query_end int,"
    "target_begin int,"
    "target_end_optimal int,"
    "evalue float,"
    "bitscore float"
    ")"
)

if pl.Path(args.output).is_file():
    os.remove(args.output)

command = ["sqlite3", args.output, create_table_cmd]

subprocess.run(command)

subprocess.run(["sqlite3", args.output, "CREATE INDEX id_idx on sequences(id)"])

subprocess.run(
    ["sqlite3", args.output, ".mode tabs", f".import {tmp_output_name} sequences"]
)


if pl.Path(f"{tmp_output_name}").is_file():
    os.remove(f"{tmp_output_name}")
