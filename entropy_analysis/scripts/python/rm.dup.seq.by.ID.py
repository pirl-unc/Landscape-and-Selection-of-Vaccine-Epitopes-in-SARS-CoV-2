#!/usr/local/bin/python3

# Remove duplicate sequences by ID
# Input file = fasta format

import sys
from Bio import SeqIO

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]

with open(out_fasta, 'a') as outFile:
    record_ids = list()
    for record in SeqIO.parse(in_fasta, 'fasta'):
        if record.id not in record_ids:
            record_ids.append( record.id)
            SeqIO.write(record, outFile, 'fasta')
