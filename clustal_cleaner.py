"""

CLUSTAL Cleaner

Clear the CLUSTAL output and position the forward and reverse primers
(c) Olivier Friard

"""

__version__ = '3'
__version_date__ = "2021-10-19"

ROW_HEADER_ID = 5

import sys
import pathlib as pl
import argparse

from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Clustal Cleaner',)

parser.add_argument('-i', action="store", dest="input")
parser.add_argument('-m', action="store", dest="model")
parser.add_argument('-s', action="store", dest="sequence_file")
parser.add_argument('-g', action="store_true", dest="group")

args = parser.parse_args()


if not args.model:
    print('Model sequence ID not specified!', file=sys.stderr)
    sys.exit(1)


if not args.input:
    print('Input file not specified!', file=sys.stderr)
    sys.exit(1)

if not pl.Path(args.input).is_file:
    print("File not found", file=sys.stderr)
    sys.exit(2)


# read sequences from alignment
align = AlignIO.read(args.input, "clustal")

ref_seq = ""
ref_id = ""
max_len_id = 0
group = {}

sequences = {}
for seq in align:

    if seq.id.upper() == args.model.upper():
        ref_seq = str(seq.seq)
        ref_id = seq.id
        max_len_id = max(max_len_id, len(seq.id))
        continue

    if args.group:
        for seq2 in sequences:
            if str(seq.seq) == sequences[seq2]:
                group[seq2] += 1
                break
        else:
            sequences[seq.id] = str(seq.seq)
            group[seq.id] = 1
            max_len_id = max(max_len_id, len(seq.id))

    else:
        sequences[seq.id] = str(seq.seq)
        max_len_id = max(max_len_id, len(seq.id))

# modify seq id for group
if args.group:
    count = 1
    for seq_id in group:
        if group[seq_id] > 1:
            sequences[f"group #{count} ({group[seq_id]} seq.)"] = sequences[seq_id]
            del sequences[seq_id]
            count += 1
    # update max_len_id
    max_len_id = max([len(id) for id in sequences])

if not ref_id:
    print(f"Model sequence ID not found in alignment!\nID list:\n{' '.join(sequences.keys())}", file=sys.stderr)
    sys.exit(1)


# verif seq len
if len(set([len(sequences[k]) for k in sequences])) != 1:
    print("error reading the CLUSTAL output. Different sequence lengths were found", file=sys.stderr)
    sys.exit(1)


# sequences to position
if args.sequence_file:
    if not pl.Path(args.sequence_file).is_file:
        print("Sequence file not found", file=sys.stderr)
        sys.exit(2)

    handle = open(args.sequence_file, "r")
    seq2position = {}
    seq_idx = {}
    polarity = {}

    for record in SeqIO.parse(handle, "fasta"):
        seq2position[record.id] = str(record.seq)
        print(f"Sequence to position {record.id} {seq2position[record.id]} ->  ", end="", file=sys.stderr)
        # search sequence
        for id in sequences:
            if seq2position[record.id] in sequences[id]:
                seq_idx[record.id] = sequences[id].index(seq2position[record.id])
                polarity[record.id] = ""
                print(f"FOUND", file=sys.stderr)
                break
        else:
            # search rev-comp
            print("NOT FOUND", file=sys.stderr)
            seq2position[record.id] = str(Seq(seq2position[record.id]).reverse_complement())
            print(f"rev-comp sequence to position {record.id} {seq2position[record.id]} ->  ", end="", file=sys.stderr)
            for id in sequences:
                if seq2position[record.id] in sequences[id]:
                    seq_idx[record.id] = sequences[id].index(seq2position[record.id])
                    polarity[record.id] = " (rev-comp)"
                    print(f"FOUND", file=sys.stderr)
                    break
            else:
                print(f"NOT FOUND", file=sys.stderr)



SPAN_F_OPEN = '<span style="background-color: yellow">'
SPAN_CLOSE = '</span>'
SPAN_R_OPEN = '<span style="background-color: lime">'

SPAN_OPEN = {"": '<span style="background-color: yellow">',
             " (rev-comp)": '<span style="background-color: lime">'}

print("<html><head></head><body><pre>")


# header
header = " " * [len(sequences[k]) for k in sequences][0]
header_out = header

# add sequences to header
if args.sequence_file:

    # pre-header for seq names
    pre_header = " " * [len(sequences[k]) for k in sequences][0]
    for seq_id in seq_idx:
        pre_header = pre_header[0:seq_idx[seq_id]] + seq_id + polarity[seq_id] + header[seq_idx[seq_id] + len(seq_id + polarity[seq_id]):]
    print(" " * (max_len_id + ROW_HEADER_ID), end="")
    print(pre_header)


    # position
    for seq_id in seq_idx:
        header_out = header_out[0:seq_idx[seq_id]] + seq2position[seq_id]  + header[seq_idx[seq_id] + len(seq2position[seq_id]):]

    # colorize
    for seq_id in seq_idx:
        header_out = header_out.replace(seq2position[seq_id], SPAN_OPEN[polarity[seq_id]] + seq2position[seq_id] + SPAN_CLOSE)


print(" " * (max_len_id + ROW_HEADER_ID), end="")
print(header_out)


# reference sequence
print(ref_id, end="")
print(" " * (max_len_id - len(ref_id) + ROW_HEADER_ID), end="")

ref_seq_out = ref_seq
if args.sequence_file:
    for seq_id in seq_idx:
        ref_seq_out = ref_seq_out.replace(seq2position[seq_id], SPAN_OPEN[polarity[seq_id]] + seq2position[seq_id] + SPAN_CLOSE)

print(ref_seq_out)


for id in sequences:
    seq = sequences[id]

    cleaned_seq = ""
    for idx, nt in enumerate(seq):
        if nt == ref_seq[idx]:
            cleaned_seq += "."
        else:
            cleaned_seq += nt

    if args.sequence_file:
        target = {}
        for seq_id in seq_idx:
            target[seq_id] = ""
            for idx, nt in enumerate(seq[seq_idx[seq_id]:seq_idx[seq_id] + len(seq2position[seq_id])]):
                if nt == seq2position[seq_id][idx]:
                    target[seq_id] += "."
                else:
                    target[seq_id] += nt

    print(id, end="")
    print(" " * (max_len_id - len(id) + ROW_HEADER_ID), end="")

    output = cleaned_seq
    if args.sequence_file:
        count = 0
        for seq_id in seq_idx:
            output = output[0:seq_idx[seq_id]] + (str(count) * len(target[seq_id])) + output[seq_idx[seq_id] + len(seq2position[seq_id]):]
            count += 1
            if count == 10:
                print("too many sequences to position (<10)", file=sys.stderr)
                sys.exit()

        count = 0
        for seq_id in seq_idx:
            output = output.replace(str(count) * len(target[seq_id]), SPAN_OPEN[polarity[seq_id]] + target[seq_id] + SPAN_CLOSE)
            count += 1

    print(output)

print("</pre></body></html>")
