"""
extract metadata from provision.json file from GISAID
"""

import sys
import pathlib as pl

if not pl.Path(sys.argv[1]).is_file:
    print("file not found")
    sys.exit()

with open(sys.argv[1], "r") as f_in:
    print(f_in.readline())
    

