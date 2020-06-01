#!/usr/bin/env python3

# In file = MSA in fasta format

from Bio import AlignIO
from collections import Counter
import numpy as np
import sys

input_mafft = sys.argv[1]

# Initialize list for length of sequence
length_list = list()

# Get alignment from fasta file
alignment = AlignIO.read(open(input_mafft), 'fasta')

# Get length of each sequence & add to list
for record in alignment:
    length_list.append(len(record.seq))

# Get unique length value from list
seq_length = next(iter(set(length_list)))

# Loop through each sequence position
for position in range(0, seq_length): #seq_length
    res = Counter(alignment[:, position].upper().replace("X", "").replace("-", "")) # find count of each AA in that position
    list_len = len(alignment[:, position].upper().replace("X", "").replace("-", "")) # find number of unique AAs
    sum = 0
    pair_list = list() # treat each AA & its count as a pair
    for k,v in res.items(): # Loop through each AA
        pair_list.append((k, v))
        proportion = v/list_len # find proportion
        sum = sum + (proportion * np.log(proportion)) # entropy calculation & then add to sum
    if sum != 0:
        sum = float(sum * -1)
    print(position, sum, sep = "\t")
