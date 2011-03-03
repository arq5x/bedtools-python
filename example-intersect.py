#!/usr/bin/env python

import unittest
import os

from bedtools import Interval, IntervalFile

def main():
    """
    Goal: find all rmsk annotations that overlap exons.
    
    Illustrates usage of IntervalFile.all_hits() 
    """
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")

    # Loop through each exon...
    for exon in exons:
        # search and retreive all rmsk features that overlap this exon
        for rmsk_hit in rmsk.all_hits(exon):
            # separate prints for clarity...
            print exon.chrom,     exon.start,     exon.end,     # print the exon feature
            print rmsk_hit.chrom, rmsk_hit.start, rmsk_hit.end, # print the rmsk feature
            print exon.o_start,   exon.o_end,     exon.o_amt    # print the start, end and amout of the overlap


if __name__ == "__main__":
    main()
