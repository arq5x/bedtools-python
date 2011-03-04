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
            print i.chrom, i.start, i.end, i.strand,
            print hit.chrom, hit.start, hit.end, hit.name, hit.strand,
            # o_start, o_end, o_amt are the start, end and amount of overlap
            print hit.o_start, hit.o_end, hit.o_amt

    # Example 2:
    #    Method: IntervalFile.any_hits()
    #    Report the BAM alignment if it has _any_ overlap with rmsk
    for al in bam:
        strand = "+"
        if al.is_reverse: strand = "-"
        i = Interval(bam.getrname(al.rname), al.pos, al.aend, strand)

        if rmsk.any_hits(i, same_strand=True, ovlp_pct=0.75):
            print i.chrom, i.start, i.end, i.strand

    # Example 3:
    #    Method: IntervalFile.count_hits()
    #    Report the _count_ of overlaps with rmsk for each BAM alignment
    for al in bam:
        strand = "+"
        if al.is_reverse: strand = "-"
        i = Interval(bam.getrname(al.rname), al.pos, al.aend, strand)

        hit_count = rmsk.count_hits(i, same_strand=True, ovlp_pct=0.75)
        print i.chrom, i.start, i.end, i.strand, hit_count

if __name__ == "__main__":
    main()