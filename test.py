#!/usr/bin/env python

import unittest
import os
import os.path as op

from bedtools import IntervalFile

def main():
    
    rmskFile = IntervalFile("testData/rmsk.hg18.chr21.bed")
    rmskFile.loadIntoMap()
    
    hits = rmskFile("chr21", 0, 100000, "+", False)

    for hit in hits():
        print hit.start

if __name__ == "__main__":
    main()