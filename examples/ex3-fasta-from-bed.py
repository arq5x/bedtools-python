#!/usr/bin/env python
import os
import string
from pysam    import Fastafile
from bedtools import IntervalFile

def main():
    # setup a reverse_complement translation
    rev_table=string.maketrans('ACGTacgt', 'TGCAtgca')
    def revcomp(seq, rev_table):
        return seq.translate(rev_table)
        
    # open your fasta file
    fasta  = Fastafile("bedtools/tests/data/chr21.fa")
    # open your bed file
    bed    = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    
    # for each bed, grab the the DNA in that interval 
    for b in bed:
        # grab the seq, rev. comp if necessary
        seq = fasta.fetch(b.chrom, b.start, b.end)  
        if b.strand == "-":
            seq = revcomp(seq, rev_table)
        # print the interval and the seq
        print b.chrom, b.start, b.end, b.strand, seq

if __name__ == "__main__":
    main()