"""
align a query sequence with the smith-waterman algorithm (skbio) with sequences in a multi fasta file
(c) Olivier Friard
"""

__version__ = '3'
__version_date__ = "2021-05-20"

import os
from multiprocessing import Pool, cpu_count
from skbio.alignment import StripedSmithWaterman
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
import argparse
import pathlib as pl

MIN_ALIGN_LENGTH = 0.5
MIN_IDENTITY = 0.5

MATCH_SCORE = 5
MISMATCH_SCORE = -4
GAP_OPEN_PENALTY = 12
GAP_EXTEND_PENALTY = 4


'''
alignment.aligned_query_sequence      alignment.query_begin                 alignment.target_begin
alignment.aligned_target_sequence     alignment.query_end                   alignment.target_end_optimal
alignment.cigar                       alignment.query_sequence              alignment.target_end_suboptimal
alignment.is_zero_based(              alignment.set_zero_based(             alignment.target_sequence
alignment.optimal_alignment_score     alignment.suboptimal_alignment_score  
'''


def output(id, description, alignment):
    return f"{id}\t{description}\t{alignment['frame']}\t{alignment['identity']}\t{alignment['score']}\t{alignment['align_length']}\t{alignment['target_length']}\t{alignment['aligned_query_sequence']}\t{alignment['aligned_target_sequence']}\t{alignment['query_begin']}\t{alignment['query_end']}\t{alignment['target_begin']}\t{alignment['target_end_optimal']}"


def align(input): # target_seq, target_id, target_description):

    target_seq, target_id, target_description = input
    result = ""

    frame = "f"

    idx = -1
    found_100 = False
    while True:
        idx = target_seq.find(query_sequence, idx + 1)
        if idx == -1:
            break
        found_100 = True
        alignment_f = {"frame": frame,
                       "score": query_sequence_length * MATCH_SCORE,
                       "identity": "100.00",
                       "align_length": len(query_sequence),
                       "target_length": len(target_seq),
                       "aligned_query_sequence": query_sequence,
                       "aligned_target_sequence": query_sequence,
                       "query_begin": 1,
                       "query_end": len(query_sequence),
                       "target_begin": idx + 1,
                       "target_end_optimal": idx  + len(query_sequence)
                      }
        result += output(target_id, target_description, alignment_f)

    if not found_100:
        ssw_alignment = query(target_seq)
        assert len(ssw_alignment.aligned_query_sequence) == len(ssw_alignment.aligned_target_sequence)

        # identity
        matches = 0
        for i, j in zip(ssw_alignment.aligned_query_sequence, ssw_alignment.aligned_target_sequence):
            if i == j:
                matches += 1
        identity = 100 * matches / len(ssw_alignment.aligned_query_sequence)

        alignment_f = {"frame": frame,
                       "score": ssw_alignment['optimal_alignment_score'],
                       "identity": round(identity, 2),
                       "align_length": len(ssw_alignment.aligned_target_sequence),
                       "target_length": len(record.seq),
                       "aligned_query_sequence": ssw_alignment.aligned_query_sequence,
                       "aligned_target_sequence": ssw_alignment.aligned_target_sequence,
                       "query_begin":  ssw_alignment.query_begin,
                       "query_end": ssw_alignment.query_end,
                       "target_begin": ssw_alignment.target_begin,
                       "target_end_optimal": ssw_alignment.target_end_optimal,
                       }
        if MIN_IDENTITY or MIN_ALIGN_LENGTH:
            if alignment_f["identity"] >= MIN_IDENTITY and alignment_f["align_length"] / query_sequence_length >= MIN_ALIGN_LENGTH:
                result += output(target_id, target_description, alignment_f)
        else:
            result += output(record.id, record.description, alignment_f)


    frame = "r"
    idx = -1
    found_100 = False
    while True:
        idx = target_seq.find(query_sequence_revcomp, idx + 1)
        if idx == -1:
            break
        found_100 = True

        alignment_r = {"frame": frame,
                       "score": query_sequence_length * MATCH_SCORE,
                       "identity": "100.00",
                       "align_length": query_sequence_length,
                       "target_length": len(target_seq),
                       "aligned_query_sequence": query_sequence,
                       "aligned_target_sequence": query_sequence,
                       "query_begin": 1,
                       "query_end": len(query_sequence_revcomp),
                       "target_begin": idx + 1,
                       "target_end_optimal": idx + len(query_sequence_revcomp)
                       }

        result += output(target_id, target_description, alignment_r)

    if not found_100:
        ssw_alignment_revcomp = query_revcomp(target_seq)
        # identity
        matches = 0
        for i, j in zip(ssw_alignment_revcomp.aligned_query_sequence, ssw_alignment_revcomp.aligned_target_sequence):
            if i == j:
                matches += 1
        identity = 100 * matches / len(ssw_alignment_revcomp.aligned_query_sequence)

        alignment_r = {"frame": frame,
                       "score": ssw_alignment_revcomp['optimal_alignment_score'],
                       "identity": round(identity, 2),
                       "align_length": len(ssw_alignment_revcomp.aligned_target_sequence),
                       "target_length": len(record.seq),
                       #"aligned_query_sequence": ssw_alignment_revcomp.aligned_query_sequence,
                       "aligned_query_sequence": str(Seq(ssw_alignment_revcomp.aligned_query_sequence).reverse_complement()),
                       #"aligned_target_sequence": ssw_alignment_revcomp.aligned_target_sequence,
                       "aligned_target_sequence": str(Seq(ssw_alignment_revcomp.aligned_target_sequence).reverse_complement()),
                       "query_begin":  ssw_alignment_revcomp.query_begin,
                       "query_end": ssw_alignment_revcomp.query_end,
                       "target_begin": ssw_alignment_revcomp.target_begin,
                       "target_end_optimal": ssw_alignment_revcomp.target_end_optimal,
                      }

        if MIN_IDENTITY or MIN_ALIGN_LENGTH:
            if alignment_r["identity"] >= MIN_IDENTITY and alignment_r["align_length"] / query_sequence_length >= MIN_ALIGN_LENGTH:
               result += output(record.id, record.description, alignment_r)
        else:
            result += output(record.id, record.description, alignment_r)

    return result

def align_mp(seq_list):

    with Pool(n_cpu) as p:
        results = p.map(align, seq_list)
    return results


def main():

    # header
    print(("id\tdescription\tframe\tidentity\tscore\talign_length\ttarget_length\t"
           "aligned_query_sequence\taligned_target_sequence\tquery_begin\tquery_end\t"
           "target_begin\ttarget_end_optimal"), file=output_file)

    seq_list = []
    count = 0
    for record in SeqIO.parse(target_file, "fasta"):
        count += 1
        seq_list.append((str(record.seq), record.id, record.description))
        if count == 10000:
            count = 0
            results = align_mp(seq_list)
            for result in results:
                print(result, file=output_file)
            seq_list = []

    # check if seq_list is not empty
    if seq_list:
        results = align_mp(seq_list)
        for result in results:
            print(result, file=output_file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog="SW",
                                    usage="\nsw.py -q QUERY_PATH -t TARGET_PATH -c N_CORES -o OUTPUT_PATH",
                                    description="Alignment with Smith-Waterman")
    parser.add_argument("-q", "--query",  action="store", dest="query", type=str, help="Path of the query file (FASTA format)")
    parser.add_argument("-t", "--target",  action='store', dest="target", type=str, help="Path of the target sequences (FASTA format)")
    parser.add_argument("-c", "--cpu", action="store", dest="cpu", default=16, type=int, help="Set number of CPU/cores to use (default all)")
    parser.add_argument("-o", "--output", action="store", dest="output", help="Set path for the output file")
    parser.add_argument("-v", "--version", action='version', version=f"%(prog)s v.{__version__} {__version_date__}")
    parser.add_argument("--min-align-len", action='store', dest="min_align_len", type=float, help="Minimal length of alignment (0-100% of query length)")
    parser.add_argument("--min-identity", action='store', dest="min_identity", type=float, help="Minimal identity (0-100%)")

    MIN_ALIGN_LENGTH = 0.5
    MIN_IDENTITY = 0.5

    args = parser.parse_args()

    if not args.query:
        print('Query file not found!')
        sys.exit(1)
    else:
        query_file = args.query

    if not args.target:
        print('Target file not found!')
        sys.exit(1)
    else:
        target_file = args.target

    if not args.cpu:
        n_cpu = cpu_count()
    else:
        n_cpu = int(args.cpu)

    if args.output:
        if pl.Path(args.output).is_file():
            os.remove(args.output)
        output_file = open(args.output, "a")
    else:
        output_file = sys.stdout

    if args.min_align_len:
        MIN_ALIGN_LENGTH = args.min_align_len / 100
    else:
        MIN_ALIGN_LENGTH = 0

    if args.min_identity:
        MIN_IDENTITY = args.min_identity / 100
    else:
        MIN_IDENTITY = 0


    # read the query sequence from argv #1
    for record in SeqIO.parse(query_file, "fasta"):
        query_id = record.id
        query_sequence = str(record.seq)
        query_sequence_length = len(query_sequence)

    query = StripedSmithWaterman(query_sequence, zero_index=True, gap_open_penalty=GAP_OPEN_PENALTY, gap_extend_penalty=GAP_EXTEND_PENALTY, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE)

    query_sequence_revcomp = str(Seq(query_sequence).reverse_complement())
    query_revcomp = StripedSmithWaterman(query_sequence_revcomp, zero_index=True, gap_open_penalty=GAP_OPEN_PENALTY, gap_extend_penalty=GAP_EXTEND_PENALTY, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE)

    main()

