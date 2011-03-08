#!/usr/bin/env python

import unittest
import os

from bedtools import Interval, IntervalFile

def main():
    """
    Examples of printing each interval in an interval file.
     - Works with BED, GTF and VCF files.
     - Can be uncompressed or GZIP compressed.
    """    
    
    # 0.1 Each interval in a BED file
    for exon in IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed"):
        print exon
        
    # 0.2 Each gene in a GTF file
    for gene in IntervalFile("bedtools/tests/data/genes.hg18.chr21.gtf"):
        print gene

    # 0.3 Each gene in a _compressed_ GTF file
    for gene in IntervalFile("bedtools/tests/data/genes.hg18.chr21.gtf.gz"):
        print gene
        



if __name__ == "__main__":
    main()
