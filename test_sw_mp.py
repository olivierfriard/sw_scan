import os
from multiprocessing import Pool
from skbio.alignment import StripedSmithWaterman
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime


query = StripedSmithWaterman("GGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGAC")

def align(seq2):
    #print(seq2)
    # query = StripedSmithWaterman(seq1)
    alignment = query(seq2)
    return str(alignment)


# print(align("acatcatcat", "acatcatcat"))

def start_align(l):
    with Pool(os.cpu_count()) as p:
        results = p.map(align, l)
    return results


def main():

    handle = open(sys.argv[1], "r")
    main_results = []
    list1 = []
    for record in SeqIO.parse(handle, "fasta"):
        list1.append(str(record.seq))

        if len(list1) > 10000:
            results = start_align(list1)
            main_results.extend(results)
            list1 = []


    print(main_results)

if __name__ == '__main__':
    main()

