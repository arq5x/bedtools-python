#!/usr/bin/env python

import unittest
import os

from bedtools import IntervalFile, Overlap

def main():

    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    
    # find exons that overlap with repeat annotations
    for ex in exons:
        # retrieve repeats that overlap this exon
        for hit in rmsk.search(ex.chrom, ex.start, ex.end):
            # report the _full_ features that overlap
            print ex.chrom,  ex.start,  ex.end,  ex.name,\
                  hit.chrom, hit.start, hit.end, hit.name,
            # extract the coordinates of the overlap, 
            # as well as the amount of overlap in b.p.
            o = Overlap(ex.start, ex.end, hit.start, hit.end)
            print "overlap_start=" + str(o.overlap_start),\
                  " overlap_end="  + str(o.overlap_end),\
                  " overlap_amt="  + str(o.overlap_amt)

if __name__ == "__main__":
    main()
