#!/usr/bin/env python

import unittest
import os

from bedtools import IntervalFile

def main():

    rmskFile = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    rmskFile.loadIntoMap()
    
    print "* either strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768)
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other

    print "* + strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, "+")
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other
        
    print "* either strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, None, 0.1)
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other
        
    print "* - strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, "-", 0.1)
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other
        
    
    rmskFile = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed.2.gz")
    rmskFile.loadIntoMap()

    print "* gzip file, either strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768)
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other

    print "* gzip file, + strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, "+")
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other

    print "* gzip file, either strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, None, 0.1)
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other

    print "* gzip file, - strand, any overlap *"
    hits = rmskFile.findOverlaps("chr21", 9719768, 9739768, "-", 0.1)
    for hit in hits:
        print hit.chrom, hit.start, hit.end, hit.strand, hit.name, hit.other

if __name__ == "__main__":
    main()
