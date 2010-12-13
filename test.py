#!/usr/bin/env python

import unittest
import os

from bedtools import IntervalFile

def main():

    rmskFile = IntervalFile("testData/rmsk.hg18.chr21.bed")
    rmskFile.loadIntoMap()

    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, "+", False)

    for hit in hits:
        print hit.start, hit.end
        # print string() stuff not ok??
        print hit.chrom, hit

if __name__ == "__main__":
    main()
