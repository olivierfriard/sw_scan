"""

extract metadata from provision.json file from GISAID

create a fasta one-row file


example of metadata:

covsurver_prot_mutations: (Spike_D614G,NSP15_A283V,NSP12_P323L)
covsurver_uniquemutlist: 
covv_accession_id: EPI_ISL_426900
covv_clade: G
covv_collection_date: 2020-03-25
covv_host: Human
covv_lineage: B.1
covv_location: Oceania / Australia / Northern Territory
covv_passage: Original
covv_subm_date: 2020-04-17
covv_type: betacoronavirus
covv_variant: 
covv_virus_name: hCoV-19/Australia/NT12/2020
gc_content: 0.3796742758876488
is_complete: True
is_high_coverage: True
is_reference: False
n_content: 0.006864911928203068
pangolin_lineages_version: 2021-09-16
sequence_length: 29862



"""

import sys
import json
import pathlib as pl

if not pl.Path(sys.argv[1]).is_file:
    print("file not found")
    sys.exit()

with open(sys.argv[1], "r") as f_in, open(sys.argv[2], "w") as f_out:
    for line in f_in:
        #print(line)
        line_from_json = json.loads(line)
        for k in list(line_from_json.keys()):

            # one-row fasta file
            if k == 'sequence':
                print(f">{line_from_json['covv_accession_id']}", file=f_out)
                print(line_from_json[k].replace("\n", ""), file=f_out)
                continue

            # metadata
            # print(f"{k}: {line_from_json[k]}")
            print(f"{line_from_json[k]}", end="\t")

        print()


