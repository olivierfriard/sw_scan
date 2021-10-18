

forward_seq = "AGCTTGCCAAGATCGATGGCGTGACCAATTTT"


import sys
import pathlib as pl

from Bio import AlignIO

if not pl.Path(sys.argv[1]).is_file:
    print("File not found")
    sys.exit(2)

align = AlignIO.read(sys.argv[1], "clustal")

ref_seq = ""
ref_id = ""
max_len_id = 0

sequences = {}

for seq in align:
    if not ref_seq:
        ref_seq = str(seq.seq)
        ref_id = seq.id
        max_len_id = len(ref_id)
        continue

    sequences[seq.id] = str(seq.seq)
    max_len_id = max(max_len_id, len(seq.id))

    # print(len(str(seq.seq)))

# verif seq len
if len(set([len(sequences[k]) for k in sequences])) != 1:
    print("error reading the CLUSTAL output. Different sequence lengths were found")
    sys.exit(1)



if forward_seq in ref_seq:
    forward_idx = ref_seq.index(forward_seq)
else:
    forward_idx = -1

SPAN_OPEN = '<span style="background-color: yellow">'
SPAN_CLOSE = '</span>'

print("<pre>")

print(ref_id, end="")
print(" " * (max_len_id - len(ref_id) + 5), end="")

if forward_idx != -1:
    ref_seq_out = ref_seq.replace(forward_seq, SPAN_OPEN + forward_seq + SPAN_CLOSE)
else:
    ref_seq_out = ref_seq

print(ref_seq_out)


for id in sequences:
    output = ""
    for idx, nt in enumerate(sequences[id]):
        if nt == ref_seq[idx]:
            output += "."
        else:
            output += nt
    print(id, end="")
    print(" " * (max_len_id - len(id) + 5), end="")

    if forward_idx != -1:
        output_out = output[0:forward_idx] + SPAN_OPEN  + output[forward_idx:forward_idx + len(forward_seq)]  + SPAN_CLOSE + output[forward_idx + len(forward_seq):]
    else:
        output_out = output

    print(output_out)

    
print("</pre>")