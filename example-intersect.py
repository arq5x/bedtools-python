#!/usr/bin/env python

import unittest
import os

from bedtools import IntervalFile

def main():

    snps = IntervalFile("bedtools/tests/data/snps.hg18.chr21.bed")
    rmsk = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    
    # find snps that overlap with repeat annotations
    for a in snps:
        for hit in rmsk.search(a.chrom, a.start, a.end):
           print a.chrom, a.start, a.end, a.name,
           print hit.chrom, hit.start, hit.end, hit.name

if __name__ == "__main__":
    main()
