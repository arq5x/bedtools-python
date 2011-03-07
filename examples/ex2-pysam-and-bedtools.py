#!/usr/bin/env python
import os

from pysam    import Samfile
from bedtools import Interval, IntervalFile

def main():

    bam  = Samfile("bedtools/tests/data/NA18152.bam", "rb")
    rmsk = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    
    # Example 1:
    #    Method: IntervalFile.all_hits()
    #    Report _all_ of the rmsk features that overlap with the BAM alignment
    for al in bam:
        strand = "+"
        if al.is_reverse: strand = "-"
        i = Interval(bam.getrname(al.rname), al.pos, al.aend, strand)
        
        for hit in rmsk.all_hits(i, same_strand=True, ovlp_pct=0.75):
            print "\t".join(str(x) for x in [i,hit])

if __name__ == "__main__":
    main()
    