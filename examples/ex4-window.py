#!/usr/bin/env python

import unittest
import os

from bedtools import Interval, IntervalFile

def main():
    """
    """

    ##########################################################
    # ex1. Report the coordinates of overlap b/w exons and rmsk
    #
    # Equivalent to: intersectBed -a exons -b rmsk 
    # Uses:           IntervalFile.all_hits()
    ##########################################################
    exons  = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk   = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    
    # allow 1kb of "slop" on each side of the exon 
    # when looking for hits
    window = 1000
    for exon in exons:
        # add the slop and search
        exon_slop = Interval(exon.chrom, exon.start-window, exon.end + window, exon.strand)
        for rmsk_hit in rmsk.all_hits(exon_slop):
            print "\t".join(str(f) for f in [exon, rmsk_hit])


if __name__ == "__main__":
    main()
