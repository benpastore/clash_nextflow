#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO

def get_seq(records, chrom, start, end) : 
    
    """ Returns genomic sequence on chrom from start to end """
    
    seq = records[chrom].seq
    
    return (str(seq)[start : end])

def get_target_sequence(df, *args) :
    
    """ Add genomic target sequence to dataframe """
    
    records = args['records']
    
    df['pad_up'] = df.apply(lambda x : x['target_start']-10 if x['target_start']-10 > 0 else 0, axis = 1 )
    
    df['pad_down'] = df.apply(lambda x : x['target_end']+10, axis = 1 )
    
    df['target_pad'] = df.apply(lambda x : get_seq(records, x['target'], x['pad_up'], x['pad_down']) , axis = 1 )
    
def prepare_rnaup(file) : 

    df = pd.read_csv(file, header = 0, sep = "\t")

    df['']