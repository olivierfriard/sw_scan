"""

CLUSTAL Cleaner

Clear the CLUSTAL output and position the forward and reverse primers
(c) Olivier Friard

"""

__version__ = '1'
__version_date__ = "2021-10-18"

import sys
import pathlib as pl
import argparse

from Bio import AlignIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(
    description='Clustal Cleaner',
)


parser.add_argument('-i', action="store", dest="input")
parser.add_argument('-m', action="store", dest="model")
parser.add_argument('-f', action="store", dest="forward")
parser.add_argument('-r', action="store", dest="reverse")
args = parser.parse_args()

forward_seq = args.forward.upper() if args.forward else ""

reverse_seq = str(Seq(args.reverse.upper()).reverse_complement()) if args.reverse else ""

if not args.input:
    print('Input file not specified!', file=sys.stderr)
    sys.exit(1)

if not pl.Path(args.input).is_file:
    print("File not found", file=sys.stderr)
    sys.exit(2)


align = AlignIO.read(args.input, "clustal")

ref_seq = ""
ref_id = ""
max_len_id = 0

sequences = {}

for seq in align:

    if seq.id.upper() == args.model.upper():
        ref_seq = str(seq.seq)
        ref_id = seq.id
        max_len_id = max(max_len_id, len(seq.id))
        continue

    sequences[seq.id] = str(seq.seq)
    max_len_id = max(max_len_id, len(seq.id))



# verif seq len
if len(set([len(sequences[k]) for k in sequences])) != 1:
    print("error reading the CLUSTAL output. Different sequence lengths were found", file=sys.stderr)
    sys.exit(1)

if forward_seq:
    for id in sequences:
        if forward_seq in sequences[id]:
            forward_idx = sequences[id].index(forward_seq)
            break
    else:
        print("forward primer not found found", file=sys.stderr)
        forward_idx = -1

if reverse_seq:
    for id in sequences:
        if reverse_seq in sequences[id]:
            reverse_idx = sequences[id].index(reverse_seq)
            break
    else:
        print("reverse primer not found found", file=sys.stderr)
        reverse_idx = -1

SPAN_F_OPEN = '<span style="background-color: yellow">'
SPAN_CLOSE = '</span>'

SPAN_R_OPEN = '<span style="background-color: lime">'


print("<pre>")


# header
header = " " * [len(sequences[k]) for k in sequences][0]

# add forward 
header_out = header[0:forward_idx] + forward_seq  + header[forward_idx + len(forward_seq):]

# add reverse
header_out = header_out[0:reverse_idx] + reverse_seq  + header_out[reverse_idx + len(reverse_seq):]

header_out = header_out.replace(forward_seq, SPAN_F_OPEN + forward_seq + SPAN_CLOSE)
header_out = header_out.replace(reverse_seq, SPAN_R_OPEN + reverse_seq + SPAN_CLOSE)

print(" " * (max_len_id + 5), end="")
print(header_out)


# reference sequence
print(ref_id, end="")
print(" " * (max_len_id - len(ref_id) + 5), end="")

if forward_idx != -1:
    ref_seq_out = ref_seq.replace(forward_seq, SPAN_F_OPEN + forward_seq + SPAN_CLOSE)
else:
    ref_seq_out = ref_seq


if reverse_idx != -1:
    ref_seq_out = ref_seq_out.replace(reverse_seq, SPAN_R_OPEN + reverse_seq + SPAN_CLOSE)


print(ref_seq_out)



for id in sequences:
    seq = sequences[id]

    cleaned_seq = ""
    for idx, nt in enumerate(seq):
        if nt == ref_seq[idx]:
            cleaned_seq += "."
        else:
            cleaned_seq += nt

    if forward_idx != -1:
        target_f = ""
        for idx, nt in enumerate(seq[forward_idx:forward_idx + len(forward_seq)]):
            if nt == forward_seq[idx]:
                target_f += "."
            else:
                target_f += nt

    if reverse_idx != -1:
        target_r = ""
        for idx, nt in enumerate(seq[reverse_idx:reverse_idx + len(reverse_seq)]):
            if nt == reverse_seq[idx]:
                target_r += "."
            else:
                target_r += nt


    print(id, end="")
    print(" " * (max_len_id - len(id) + 5), end="")

    if forward_idx != -1:
        output = cleaned_seq[0:forward_idx] + ("F" * len(target_f)) + cleaned_seq[forward_idx + len(forward_seq):]
    else:
        output = cleaned_seq

    if reverse_idx != -1:
        output = output[0:reverse_idx] + ("R" * len(target_r)) + output[reverse_idx + len(reverse_seq):]

    if forward_idx != -1:
        output = output.replace("F" * len(target_f), SPAN_F_OPEN + target_f + SPAN_CLOSE)
    if reverse_idx != -1:
        output = output.replace("R" * len(target_r), SPAN_R_OPEN + target_r + SPAN_CLOSE)

    print(output)

print("</pre>")