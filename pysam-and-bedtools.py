#!/usr/bin/env python
import os

from pysam    import Samfile
from bedtools import IntervalFile

def main():

    bam  = Samfile("bedtools/tests/data/NA18152.bam", "rb")
    rmsk = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    
    for al in bam:
        chrom = bam.getrname(al.rname)
        start = al.pos
        end   = al.aend
        name  = al.qname 
        for hit in rmsk.search(chrom, start, end):
            print chrom, start, end, name,
            print hit.chrom, hit.start, hit.end, hit.name

if __name__ == "__main__":
    main()