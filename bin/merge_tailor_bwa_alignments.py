#!/usr/bin/env python3

import os
import pandas as pd 
import argparse

def merge_bed(bwa_bed, tailor_bed, output) : 



    bwa = pd.read_csv(bwa_bed, header = None, names = ['gene_id', 'gene_start', 'gene_end', 'clash_read', 'score', 'strand', 'gene_mapped', 'smRNA_seq'], sep = "\t")
    tailor = pd.read_csv(tailor_bed, header = None, names = ['smRNA_id', 'smRNA_start', 'smRNA_end', 'clash_read', 'N_map', 'strand', 'smRNA_seq', 'tail'], sep = "\t")
    merge_df = bwa.merge(tailor, how = 'inner', on = ['clash_read', 'strand', 'smRNA_seq'])
    merge_df.to_csv(f"{output}", sep = "\t", header = True, index = False)

def get_args() : 

    """Parse command line parameters from input"""

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-bwa", type = str, required = True)
    parser.add_argument("-tailor", type = str, required = True)
    parser.add_argument("-output", type = str, required = True)

    return parser.parse_args()

def main() : 

    args = get_args()

    merge_bed(args.bwa, args.tailor, args.output)

if __name__ == "__main__" : 

    main()