#!/usr/bin/env python3

# Convert mafft alignment fasta into tab delimited file - 1st column = seqeunce ID, all other columns equal position. 
# Input = mafft alignment file

import sys

mafft_file = sys.argv[1]


with open(mafft_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            row_name = line.strip("\n").strip(">")
            print("\n", row_name, end="\t")
        else:
            line = line.strip("\n")
            line = list(line)
            for a in line:
                print(a, end="\t")
