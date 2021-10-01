"""
extract metadata from provision.json file from GISAID
"""

import sys
import json
import pathlib as pl

if not pl.Path(sys.argv[1]).is_file:
    print("file not found")
    sys.exit()

with open(sys.argv[1], "r") as f_in:
    for line in f_in:
        print(line)
        line_from_json = json.loads(line)
        print(sorted(list(line_from_json.keys())))
        break
    

