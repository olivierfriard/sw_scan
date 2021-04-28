"""
align a query sequence with the smith-waterman algorithm (skbio) with sequences in a multi fasta file
(c) Olivier Friard
"""

from skbio.alignment import StripedSmithWaterman
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

MIN_ALIGN_LENGTH = 0.5
MIN_IDENTITY = 0.5

MATCH_SCORE = 5
MISMATCH_SCORE = -4
GAP_OPEN_PENALTY = 12
GAP_EXTEND_PENALTY = 4

if len(sys.argv) < 3:
    print("no input", file=sys.stderr)
    sys.exit(1)

'''
alignment.aligned_query_sequence      alignment.query_begin                 alignment.target_begin
alignment.aligned_target_sequence     alignment.query_end                   alignment.target_end_optimal
alignment.cigar                       alignment.query_sequence              alignment.target_end_suboptimal
alignment.is_zero_based(              alignment.set_zero_based(             alignment.target_sequence
alignment.optimal_alignment_score     alignment.suboptimal_alignment_score  
'''


def output(id, description, alignment):
    print(f"{id}\t{description}\t{alignment['frame']}\t{alignment['identity']}\t{alignment['score']}\t{alignment['align_length']}\t{alignment['target_length']}\t{alignment['aligned_query_sequence']}\t{alignment['aligned_target_sequence']}\t{alignment['query_begin']}\t{alignment['query_end']}\t{alignment['target_begin']}\t{alignment['target_end_optimal']}")


# read the query sequence from argv #1
for record in SeqIO.parse(open(sys.argv[1]), "fasta"):
    query_id = record.id
    query_sequence = str(record.seq)
    query_sequence_length = len(query_sequence)

query = StripedSmithWaterman(query_sequence, zero_index=True, gap_open_penalty=GAP_OPEN_PENALTY, gap_extend_penalty=GAP_EXTEND_PENALTY, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE)

query_sequence_revcomp = str(Seq(query_sequence).reverse_complement())
query_revcomp = StripedSmithWaterman(query_sequence_revcomp, zero_index=True, gap_open_penalty=GAP_OPEN_PENALTY, gap_extend_penalty=GAP_EXTEND_PENALTY, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE)

# print(query_sequence)
# print(query_sequence_revcomp)
# sys.exit()

count = 0

print(("id\tdescription\tframe\tidentity\tscore\talign_length\ttarget_length\taligned_query_sequence\taligned_target_sequence\tquery_begin\tquery_end\ttarget_begin\ttarget_end_optimal"))

for record in SeqIO.parse(sys.argv[2], "fasta"):
    count += 1

    frame = "f"
    #if query_sequence in str(record.seq):
    target_seq = str(record.seq)
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
                       "target_length": len(record.seq),
                       "aligned_query_sequence": query_sequence,
                       "aligned_target_sequence": query_sequence,
                       "query_begin": 1,
                       "query_end": len(query_sequence),
                       "target_begin": idx + 1,
                       "target_end_optimal": idx  + len(query_sequence)
                      }
        output(record.id, record.description, alignment_f)

    if not found_100:
        ssw_alignment = query(str(record.seq))
        assert len(ssw_alignment.aligned_query_sequence) == len(ssw_alignment.aligned_target_sequence)

        # identity
        matches = 0
        for i,j in zip(ssw_alignment.aligned_query_sequence, ssw_alignment.aligned_target_sequence):
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

        if alignment_f["identity"] >= MIN_IDENTITY and alignment_f["align_length"] / query_sequence_length >= MIN_ALIGN_LENGTH:
            output(record.id, record.description, alignment_f)

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
                       "target_length": len(record.seq),
                       "aligned_query_sequence": query_sequence,
                       "aligned_target_sequence": query_sequence,
                       "query_begin": 1,
                       "query_end": len(query_sequence_revcomp),
                       "target_begin": idx + 1,
                       "target_end_optimal": idx + len(query_sequence_revcomp)
                       }

        output(record.id, record.description, alignment_r)

    if not found_100:
        ssw_alignment_revcomp = query_revcomp(str(record.seq))
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

        if alignment_r["identity"] >= MIN_IDENTITY and alignment_r["align_length"] / query_sequence_length >= MIN_ALIGN_LENGTH:
            output(record.id, record.description, alignment_r)
