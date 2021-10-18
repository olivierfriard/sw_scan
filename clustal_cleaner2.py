"""

CLUSTAL Cleaner

Clear the CLUSTAL output in correspondence of forward and reverse primers
(c) Olivier Friard

"""



import sys
import pathlib as pl
import argparse

from Bio import AlignIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser(
    description='Clustal Cleaner',
)
parser.add_argument('-i', action="store", dest="input")
parser.add_argument('-f', action="store", dest="forward")
parser.add_argument('-r', action="store", dest="reverse")
arg = parser.parse_args()

if not pl.Path(arg.input).is_file:
    print("File not found")
    sys.exit(2)

forward_seq = arg.forward.upper()


reverse_seq = str(Seq(arg.reverse.upper()).reverse_complement())
#reverse_seq = arg.reverse.upper()


align = AlignIO.read(arg.input, "clustal")

ref_seq = ""
ref_id = ""
max_len_id = 0

sequences = {}

for seq in align:

    sequences[seq.id] = str(seq.seq)
    max_len_id = max(max_len_id, len(seq.id))

    # print(len(str(seq.seq)))

# verif seq len
if len(set([len(sequences[k]) for k in sequences])) != 1:
    print("error reading the CLUSTAL output. Different sequence lengths were found")
    sys.exit(1)


for id in sequences:
    if forward_seq in sequences[id]:
        forward_idx = sequences[id].index(forward_seq)
        break
else:
    forward_idx = -1


for id in sequences:
    if reverse_seq in sequences[id]:
        reverse_idx = sequences[id].index(reverse_seq)
        break
else:
    reverse_idx = -1

SPAN_F_OPEN = '<span style="background-color: yellow">'
SPAN_CLOSE = '</span>'

SPAN_R_OPEN = '<span style="background-color: lime">'


print("<pre>")



header = " " * [len(sequences[k]) for k in sequences][0]

# add forward 
header_out = header[0:forward_idx] + forward_seq  + header[forward_idx + len(forward_seq):]

# add reverse
header_out = header_out[0:reverse_idx] + reverse_seq  + header_out[reverse_idx + len(reverse_seq):]


header_out = header_out.replace(forward_seq, SPAN_F_OPEN + forward_seq + SPAN_CLOSE)
header_out = header_out.replace(reverse_seq, SPAN_R_OPEN + reverse_seq + SPAN_CLOSE)

print(" " * (max_len_id + 5), end="")
print(header_out)


for id in sequences:
    seq = sequences[id]

    target_f = ""
    for idx, nt in enumerate(seq[forward_idx:forward_idx + len(forward_seq)]):
        if nt == forward_seq[idx]:
            target_f += "."
        else:
            target_f += nt

    target_r = ""
    for idx, nt in enumerate(seq[reverse_idx:reverse_idx + len(reverse_seq)]):
        if nt == reverse_seq[idx]:
            target_r += "."
        else:
            target_r += nt


    print(id, end="")
    print(" " * (max_len_id - len(id) + 5), end="")

    if forward_idx != -1:
        output = seq[0:forward_idx] + ("F" * len(target_f)) + seq[forward_idx + len(forward_seq):]
    else:
        output = seq

    if reverse_idx != -1:
        output = output[0:reverse_idx] + ("R" * len(target_r)) + output[reverse_idx + len(reverse_seq):]

    output = output.replace("F" * len(target_f), SPAN_F_OPEN + target_f + SPAN_CLOSE)
    output = output.replace("R" * len(target_r), SPAN_R_OPEN + target_r + SPAN_CLOSE)

    print(output)

print("</pre>")