#!/usr/bin/env python3

# Convert tab delimited file into fasta file
# Input file assumes first column is sequence ID. Every other column is postion.

import sys

in_tab_file = sys.argv[1]

with open(in_tab_file, "r") as f:
    next(f)
    for line in f:
        line_arr = line.strip("\n").split("\t")
        seq_id = ">" + line_arr[0].strip()
        del line_arr[0]
        print(seq_id)
        print("".join(line_arr))
