#!/usr/bin/python3

"""
Multi-Alignment Cleaner

(c) Olivier Friard 2021-2023

Clean multi-alignment (CLUSTAL or MUSCLE output) and position sequences from a FASTA file

-g               group the sequences that are identical
--consensus n    print a consensus sequence at n%
--stars          print a star when all nucleotides are identical
-s primers.seq   position and highlight primer/probe sequences


"""

import sys
import pathlib as pl
import argparse
import itertools
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO


__version__ = "10"
__version_date__ = "2024-05-16"

ROW_HEADER_ID = 10


IUPAC_nt_codes = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "AG": "R",
    "CT": "Y",
    "CG": "S",
    "AT": "W",
    "GT": "K",
    "AC": "M",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N",
}


IUPAC_degenerated = {
    "-": "-",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


HTML_HEADER = "<html><head></head><body><pre>"
HTML_FOOTER = "</pre></body></html>"

SPAN_OPEN = {"": '<span style="background-color: yellow">', " (rev-comp)": '<span style="background-color: lime">'}
SPAN_CLOSE = "</span>"

YELLOW = '<span style="background-color: yellow">'
GREEN = '<span style="background-color: lime">'


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Multi-Alignment Cleaner",
    )

    parser.add_argument("-v", action="store_true", dest="version", help="Print version of Multi-Alignment Cleaner")
    parser.add_argument("-i", action="store", dest="input", help="Multi-alignment file (.aln for CLUSTAL)")
    parser.add_argument("-m", action="store", dest="model", help="Reference sequence")
    parser.add_argument(
        "-s",
        action="store",
        dest="sequence_file",
        help="FASTA file containing the sequences to position on the alignment",
    )
    parser.add_argument("-g", action="store_true", dest="group", help="Group sequences when identical")
    parser.add_argument("--stars", action="store_true", dest="stars", help="Print stars when nucleotides are identical on column")
    parser.add_argument("--consensus", action="store", dest="consensus", help="Print a consensus sequence (in percent)")

    args = parser.parse_args()

    if args.version:
        print("Multi-Alignment Cleaner\n(c) Olivier Friard 2021-2023")
        print(f"v. {__version__} ({__version_date__})\n")
        sys.exit()

    if not args.model:
        print("Model sequence ID not specified!", file=sys.stderr)
        sys.exit(1)

    if not args.input:
        print("Input file not specified!", file=sys.stderr)
        sys.exit(1)

    if not pl.Path(args.input).is_file:
        print("File not found", file=sys.stderr)
        sys.exit(2)

    return args


def read_seq_from_clustal(args):
    """
    read sequences from Multi-Alignment Alignment whith BioPython
    AlignIO.read function for CLUSTAL alignment
    SeqIO.parse for aligned FASTA files
    """

    if pl.Path(args.input).suffix == ".aln":
        align = AlignIO.read(args.input, "clustal")
    else:
        align = SeqIO.parse(args.input, "fasta")

    ref_seq: str = ""
    ref_id: str = ""
    max_len_id: str = 0
    group: dict = {}
    sequences: dict = {}

    for seq in align:
        if seq.id.upper() == args.model.upper():
            ref_seq = str(seq.seq).upper()
            ref_id = seq.id
            max_len_id = max(max_len_id, len(seq.id))
            continue

        if args.group:
            for seq2 in sequences:
                if str(seq.seq) == sequences[seq2]:
                    group[seq2] += 1
                    break
            else:
                sequences[seq.id] = str(seq.seq).upper()
                group[seq.id] = 1
                max_len_id = max(max_len_id, len(seq.id))

        else:
            sequences[seq.id] = str(seq.seq).upper()
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

    error_msg = ""

    if not ref_id:
        error_msg = f"Model sequence ID not found in alignment!\nID list:\n{' '.join(sequences.keys())}"

    # verif seq len
    if len(set([len(sequences[k]) for k in sequences])) != 1:
        error_msg = "error reading the multi-alignment output. Different sequence lengths were found"

    return error_msg, sequences, ref_seq, ref_id, max_len_id


def align_sub_sequences(args, sequences):
    if not pl.Path(args.sequence_file).is_file():
        print(f"Sequence file not found: {args.sequence_file}\n", file=sys.stderr)
        sys.exit(2)

    with open(args.sequence_file, "r") as handle:
        seq2position: dict = {}
        seq_idx: dict = {}
        polarity: dict = {}

        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.description
            seq2position[seq_id] = str(record.seq)
            print(f"Sequence to position {seq_id} {seq2position[seq_id]} ->  ", end="", file=sys.stderr)
            # search sequence
            for id in sequences:
                if seq2position[seq_id] in sequences[id]:
                    seq_idx[seq_id] = sequences[id].index(seq2position[seq_id])
                    polarity[seq_id] = ""
                    print("FOUND", file=sys.stderr)
                    break
            else:
                # search rev-comp
                print("NOT FOUND", file=sys.stderr)
                seq2position[seq_id] = str(Seq(seq2position[seq_id]).reverse_complement())
                print(f"rev-comp sequence to position {seq_id} {seq2position[seq_id]} ->  ", end="", file=sys.stderr)
                for id in sequences:
                    if seq2position[seq_id] in sequences[id]:
                        seq_idx[seq_id] = sequences[id].index(seq2position[seq_id])
                        polarity[seq_id] = " (rev-comp)"
                        print("FOUND", file=sys.stderr)
                        break
                else:
                    print("NOT FOUND", file=sys.stderr)

    return seq2position, seq_idx, polarity


def reference_seq(args, ref_seq, ref_id, max_len_id, seq2position, seq_idx, polarity):
    # reference sequence
    o = ref_id + (" " * (max_len_id - len(ref_id) + ROW_HEADER_ID))

    stars = list(ref_seq)
    if args.consensus:
        # initialize consensus with reference sequence
        consensus = {}
        for idx, nt in enumerate(ref_seq):
            if idx not in consensus:
                consensus[idx] = {}

            # dispatch nt from degenerated code
            for nt2 in IUPAC_degenerated[nt]:
                if nt2 not in consensus[idx]:
                    consensus[idx][nt2] = 1 / len(IUPAC_degenerated[nt])
                else:
                    consensus[idx][nt2] += 1 / len(IUPAC_degenerated[nt])

    else:
        consensus = {}

    ref_seq_out = ref_seq
    if args.sequence_file:
        for idx, seq_id in enumerate(seq_idx):
            if idx == 0 or idx == 2:
                color = YELLOW
            elif idx == 1:
                color = GREEN
            else:
                color = YELLOW

            ref_seq_out = ref_seq_out.replace(seq2position[seq_id], color + seq2position[seq_id] + SPAN_CLOSE)
            # ref_seq_out = ref_seq_out.replace(seq2position[seq_id], SPAN_OPEN[polarity[seq_id]] + seq2position[seq_id] + SPAN_CLOSE)

    o += ref_seq_out

    return o, consensus, stars


def display_sub_sequences(sequences, seq2position, seq_idx, polarity, max_len_id):
    # header
    header = " " * [len(sequences[k]) for k in sequences][0]
    header_out = header

    o = ""
    # pre-header for seq names
    pre_header = " " * [len(sequences[k]) for k in sequences][0]
    for seq_id in seq_idx:
        pre_header = (
            pre_header[0 : seq_idx[seq_id]] + seq_id + polarity[seq_id] + header[seq_idx[seq_id] + len(seq_id + polarity[seq_id]) :]
        )
    # print(" " * (max_len_id + ROW_HEADER_ID), end="")
    o += " " * (max_len_id + ROW_HEADER_ID)

    # print(pre_header)
    o += pre_header

    # position
    for seq_id in seq_idx:
        header_out = header_out[0 : seq_idx[seq_id]] + seq2position[seq_id] + header[seq_idx[seq_id] + len(seq2position[seq_id]) :]

    # colorize
    for idx, seq_id in enumerate(seq_idx):
        if idx == 0 or idx == 2:
            color = YELLOW
        elif idx == 1:
            color = GREEN
        else:
            color = YELLOW

        header_out = header_out.replace(seq2position[seq_id], color + seq2position[seq_id] + SPAN_CLOSE)
        # header_out = header_out.replace(seq2position[seq_id], SPAN_OPEN[polarity[seq_id]] + seq2position[seq_id] + SPAN_CLOSE)

    o += "\n" + (" " * (max_len_id + ROW_HEADER_ID)) + header_out

    return o


def clean_sequences(args, seq_idx, sequences, seq2position, ref_seq, consensus, stars, max_len_id, polarity):
    o = ""

    for id in sequences:
        seq = sequences[id]

        cleaned_seq = ""
        for idx, nt in enumerate(seq):
            if nt in IUPAC_degenerated[ref_seq[idx]] and nt != "-":
                cleaned_seq += "."
            else:
                cleaned_seq += nt

            # *
            if args.stars:
                if nt != stars[idx]:
                    stars[idx] = " "

            # consensus
            if args.consensus:
                # dispatch nt from degenerated code
                for nt2 in IUPAC_degenerated[nt]:
                    if nt2 not in consensus[idx]:
                        consensus[idx][nt2] = 1 / len(IUPAC_degenerated[nt])
                    else:
                        consensus[idx][nt2] += 1 / len(IUPAC_degenerated[nt])

        if args.sequence_file:
            target = {}
            for seq_id in seq_idx:
                target[seq_id] = ""
                for idx, nt in enumerate(seq[seq_idx[seq_id] : seq_idx[seq_id] + len(seq2position[seq_id])]):
                    # if nt == seq2position[seq_id][idx]:
                    if (nt != "-") and (nt in IUPAC_degenerated[seq2position[seq_id][idx]]):
                        target[seq_id] += "."
                    else:
                        target[seq_id] += nt

        # print(id, end="")
        o += id
        # print(" " * (max_len_id - len(id) + ROW_HEADER_ID), end="")

        o += " " * (max_len_id - len(id) + ROW_HEADER_ID)

        output = cleaned_seq
        if args.sequence_file:
            count = 0
            for seq_id in seq_idx:
                output = (
                    output[0 : seq_idx[seq_id]] + (str(count) * len(target[seq_id])) + output[seq_idx[seq_id] + len(seq2position[seq_id]) :]
                )
                count += 1
                if count == 10:
                    print("too many sequences to position (<10)", file=sys.stderr)
                    sys.exit()

            count = 0
            for idx, seq_id in enumerate(seq_idx):
                if idx == 0 or idx == 2:
                    color = YELLOW
                elif idx == 1:
                    color = GREEN
                else:
                    color = YELLOW

                output = output.replace(str(count) * len(target[seq_id]), color + target[seq_id] + SPAN_CLOSE)
                # output = output.replace(str(count) * len(target[seq_id]), SPAN_OPEN[polarity[seq_id]] + target[seq_id] + SPAN_CLOSE)
                count += 1

        # print(output)
        o += output + "\n"

    return o


def make_consensus(args, consensus, max_len_id: int) -> str:
    o = ""
    row_header = f"{args.consensus}%"

    o += row_header

    o += " " * (max_len_id - len(row_header) + ROW_HEADER_ID)

    for idx in sorted(consensus.keys()):
        print(f"{idx=}", file=sys.stderr)

        print(f"{consensus[idx]=}", file=sys.stderr)

        total = sum([consensus[idx][k] for k in consensus[idx]])

        for nt in consensus[idx]:
            if consensus[idx][nt] / total >= int(args.consensus) / 100:
                o += nt
                break
        else:
            nt_list = list(consensus[idx])
            print(f"{nt_list=}", file=sys.stderr)

            # remove gap '-'
            # nt_list.remove('-')

            flag_stop = False
            # test
            for n in range(1, 4 + 1):
                print(f"{n=}", file=sys.stderr)

                for c in itertools.combinations(nt_list, n):
                    print(f"{c=}", file=sys.stderr)

                    if sum([consensus[idx][k] for k in c]) / total >= int(args.consensus) / 100:
                        if "-" in c:
                            c_list = list(c)
                            c_list.remove("-")

                            print(f"{c_list=}", file=sys.stderr)
                            consensus_nt = IUPAC_nt_codes["".join(sorted(c_list))]
                            o += f"<u>{consensus_nt}</u>"
                        else:
                            o += IUPAC_nt_codes["".join(sorted(c))]

                        flag_stop = True
                        break

                if flag_stop:
                    break

    return o


def display_stars(args, stars, max_len_id: int) -> str:
    o = ""
    # print *
    if args.stars:
        o += " " * (max_len_id + ROW_HEADER_ID)
        for nt in stars:
            if nt == " ":
                o += nt
            else:
                o += "*"
    return o


def display(HTML_HEADER, CONSENSUS, PRIMERS_PROBES, REF_SEQ, SEQUENCES, STARS, HTML_FOOTER):
    print(HTML_HEADER)

    print(f"Multi-Alignment Cleaner v. {__version__} ({__version_date__}) by Olivier Friard<br><br>")

    print(CONSENSUS)

    print(PRIMERS_PROBES)

    print(REF_SEQ)

    print(SEQUENCES)

    print(STARS)

    print(HTML_FOOTER)


def main():
    args = parse_arguments()

    error_msg, sequences, ref_seq, ref_id, max_len_id = read_seq_from_clustal(args)

    if error_msg:
        print(error_msg, file=sys.stderr)
        sys.exit()

    # sequences to position
    if args.sequence_file:
        seq2position, seq_idx, polarity = align_sub_sequences(args, sequences)
        PRIMERS_PROBES = display_sub_sequences(sequences, seq2position, seq_idx, polarity, max_len_id)
    else:
        seq2position, seq_idx, polarity = None, None, None
        PRIMERS_PROBES = ""

    REF_SEQ, consensus, stars = reference_seq(args, ref_seq, ref_id, max_len_id, seq2position, seq_idx, polarity)

    SEQUENCES = clean_sequences(args, seq_idx, sequences, seq2position, ref_seq, consensus, stars, max_len_id, polarity)

    if args.stars:
        STARS = display_stars(args, stars, max_len_id)
    else:
        STARS = ""

    if args.consensus:
        CONSENSUS = make_consensus(args, consensus, max_len_id)
    else:
        CONSENSUS = ""

    display(HTML_HEADER, CONSENSUS, PRIMERS_PROBES, REF_SEQ, SEQUENCES, STARS, HTML_FOOTER)


if __name__ == "__main__":
    main()
