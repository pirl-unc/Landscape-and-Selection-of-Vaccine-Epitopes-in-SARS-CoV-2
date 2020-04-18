#!/usr/bin/env python

import argparse
import re

def parse_header(header):
    """
    Indexes tab-delimited string and returns dictionary with mapping of fields of interest.
    Ideally, the resulting dictionary is used to determine which column within
    a row contains a value of interest. For example, instead of assuming the
    5th position of each row contains the "Pos" value (your_row[5]), one can
    instead use your_row[hdr_idx['Pos']].

    Args:
        header (str): Tab-delimited string.

    Returns:
        dict where keys are columns of interest and values are their 0-th index in the header.
    """
    header = header.strip().split('\t')
    hdr_idx = {}
    for i in ['Pos', 'Peptide', 'ID', 'Min_nM', 'Frequency']:
        hdr_idx[i] = header.index(i)
    return hdr_idx


def det_lo_entropy(entropies, start, stop, mode, threshold):
    """
    Determines which epitopes are low or high entropy.

    Args:
        entropies (dict): Position-indexed entropies
        start (int): Start position of  epitope
        stop (int): Stop position of  epitope (inclusive)
        mode (str): Either 'avg' or 'all'. Determines what epitope-specific
                    value is compared to threshold.
        threshold (float): Threshold compared against to determine if eipotope
                           is low or high entropy.

    Returns:
        1 if low entropy, otherwise 0.
    """
    lo_entropy = 'foo'
    print(start, stop)
    print(range(start, stop))
    entrs = [float(entropies[i]) for i in range(start, stop)]
    #print(entrs)
    if mode == 'avg':
        if (sum(entrs) + 0.0)/len(entrs) < float(threshold):
            lo_entropy = '1'
        else:
            lo_entropy = '0'
    elif mode == 'all':
        if all([float(i) < float(threshold) for i in entrs]):
            lo_entropy = '1'
        else:
            lo_entropy = '0'
    return lo_entropy


def main():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--class-1-input', help="MHC Class I input file", required=True)
    parser.add_argument('-2', '--class-2-input', help="MHC Class II input file", required=True)
    parser.add_argument('-f', '--entropy-file', help="Entropy file")
    parser.add_argument('-t', '--entropy-threshold', default=0.1, help="Entropy threshold")
    parser.add_argument('-m', '--entropy-mode', default='all', help="Entropy mode (avg or all)")
    parser.add_argument('-o', '--output', help="Output file")

    args = parser.parse_args()

    # Creating header for output file.
    newf = []
    newf.append(['gene', 'c1_peptide', 'c2_peptide', 'c1_start', 'c2_start', 'c1_tot_freq', 'c2_tot_freq', 'c1_min_nm', 'c2_min_nm', 'lo_entropy'])


    # Gene start positions based off of counting AAs within MT072688.1
    prot_starts = {'orf1ab': 1,
                   'S': 7097,
                   'ORF3a': 8370,
                   'E': 8645,
                   'M': 8720,
                   'ORF6': 8942,
                   'ORF7a': 9003,
                   'ORF8': 9124,
                   'N': 9245,
                   'ORF10': 9664}

    # Initializing class II dictionary.
    class_2 = {}

    with open(args.class_2_input) as ifo:
        # Create header index
        hdr_idx =  parse_header(ifo.readline())
        print(hdr_idx)
        for line in ifo.readlines():
            line = line.strip().split('\t')
            id = line[hdr_idx['ID']].split('_')[0]
  
            # If identifier is 'orflab' then change to correct 'orf1ab'
            if id == 'orflab':
                id = 'orf1ab'

            # Initializing identifier in class II dictionary if not already present.
            if id not in class_2.keys():
                class_2[id] = {}

            # Extract information of interest from each line.
            peptide = line[hdr_idx['Peptide']]
            freq = line[hdr_idx['Frequency']]
            nm = line[hdr_idx['Min_nM']]
            start = line[hdr_idx['Pos']]
            class_2[id][peptide] = {}
            class_2[id][peptide]['freq'] = freq
            class_2[id][peptide]['nm'] = nm
            class_2[id][peptide]['start'] = start


    # Parse and extract entropy values. Note that entropies are 1-indexed and
    # assumed to only cover intragenic codons.
    entropies = {}
    with open(args.entropy_file) as efo:
        efo.readline() #Dump header for now
        for line_idx, line in enumerate(efo.readlines()):
            line = line.strip().split(',')
            entropies[line_idx+1] = line[1]

    print(sorted(entropies.keys())[:10])
    print(sorted(entropies.keys())[-1])




    with open(args.class_1_input) as ifo:
        hdr_idx =  parse_header(ifo.readline())

        for line in ifo.readlines():
            line = line.strip().split('\t')
            gen_line = []
            id = line[hdr_idx['ID']].split('_')[0]
   
            # If identifier is 'orflab' then change to correct 'orf1ab'
            if id == 'orflab':
                id = 'orf1ab'

            # Seeding with unlikely values for a peptide in order to find the best option.
            # local_c2 is a bad name here, but basically, it's the best pairing
            # for this class I epitope.
            local_c2 = {'peptide': '',
                        'start': 0,
                        'nm': 1000,
                        'freq': 0}
            if id in class_2.keys():
                for c2e in class_2[id].keys():
                    if re.search(line[hdr_idx['Peptide']], c2e):
                        # Find the highest frequency overlapping class II epitope.
                        if class_2[id][c2e]['freq'] > local_c2['freq']:
                            local_c2['peptide'] = c2e
                            local_c2['start'] = class_2[id][c2e]['start']
                            local_c2['nm'] = class_2[id][c2e]['nm']
                            local_c2['freq'] = class_2[id][c2e]['freq']

            #Checking entropy values here.
            # Renaming everything to make variable usage easier
            c1_peptide = line[hdr_idx['Peptide']]
            c1_abs_start = prot_starts[id] + int(line[hdr_idx['Pos']])
            c1_len = len(c1_peptide)
            c1_abs_stop = c1_abs_start + c1_len
            c1_freq = line[hdr_idx['Frequency']]
            c1_nm = line[hdr_idx['Min_nM']]

            c2_peptide = local_c2['peptide']
            # Need  '- 1' on following line to put NetMHCpan and NetMHCIIpan coordinates on same index.
            c2_abs_start = prot_starts[id] + int(local_c2['start']) - 1 
            c2_len = len(c2_peptide)
            c2_abs_stop = c2_abs_start + c2_len
            c2_freq = local_c2['freq']
            c2_nm = local_c2['nm']

            # Entropy is determiend using the class II epitope since class I is
            # a substring of it. Basically, if the class II epitope is low
            # entropy, then so is the class I epitope.
            if local_c2['peptide']:
                lo_entropy = ''
                print(id)
                print(c2_peptide)
                if c1_len > c2_len or c1_len == c2_len:
                    lo_entropy = det_lo_entropy(entropies, c1_abs_start, c2_abs_stop, args.entropy_mode, args.entropy_threshold)
                elif c2_len > c1_len:
                    lo_entropy = det_lo_entropy(entropies, c2_abs_start, c2_abs_stop, args.entropy_mode, args.entropy_threshold)

                # Making output line
                gen_line.append(id)
                gen_line.append(c1_peptide)
                gen_line.append(c2_peptide)
                gen_line.append(c1_abs_start)#prot_starts[id] + int(line[hdr_idx['Pos']]))
                gen_line.append(c2_abs_start)#prot_starts[id] + int(local_c2['start']))
                gen_line.append(c1_freq)#line[hdr_idx['Frequency']])
                gen_line.append(c2_freq)#local_c2['freq'])
                gen_line.append(c1_nm)#line[hdr_idx['Min_nM']])
                gen_line.append(c2_nm)#local_c2['nm'])
                gen_line.append(lo_entropy)
                #print(gen_line)
                newf.append(gen_line)

    # Write all newly generated output file lines to output file.
    with open(args.output, 'w') as ofo:
        for line in newf:
            line = [str(i) for i in line]
            ofo.write(','.join(line) + '\n')

main()
