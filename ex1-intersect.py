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
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        for rmsk_hit in rmsk.all_hits(exon):
            print "\t".join(str(f) for f in [exon.chrom, rmsk_hit.o_start, rmsk_hit.o_end])



    ##########################################################
    # ex2. Report the original features for overlapping 
    #    exons and rmsk
    #
    # Equivalent to: intersectBed -a exons -b rmsk -wa -wb
    # Uses:           IntervalFile.all_hits()
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        for rmsk_hit in rmsk.all_hits(exon):
            print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, 
                                             exon.name, exon.score, exon.strand,
                                             rmsk_hit.chrom, rmsk_hit.start, rmsk_hit.end, 
                                             rmsk_hit.name, rmsk_hit.score, rmsk_hit.strand])



    ##########################################################
    # ex3. Report the count of rmsk overlapping each exon
    #
    # Equivalent to: intersectBed -a exons -b rmsk -c
    # Uses:           IntervalFile.count_hits()
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        # get the number of hits in rmsk
        num_hits = rmsk.count_hits(exon)  
        print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, 
                                         exon.name, exon.score, exon.strand,
                                         num_hits])



    ##########################################################
    # ex4. Report exons that overlap at least one rmsk
    #
    # Equivalent to: intersectBed -a exons -b rmsk -u
    # Uses:           IntervalFile.any_hits()
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        # does this exon overlap any rmsk?
        if rmsk.any_hits(exon):
            print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, 
                                             exon.name, exon.score, exon.strand])



    ##########################################################
    # ex5. Report exons that DO NOT overlap at least one rmsk
    #
    # Equivalent to: intersectBed -a exons -b rmsk -v
    # Uses:           IntervalFile.any_hits()
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        # does this exon overlap any rmsk?
        if not rmsk.any_hits(exon):
            print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, 
                                             exon.name, exon.score, exon.strand])



    ##########################################################
    # ex6. Report overlap b/w exons and rmsk on the same strand
    #
    # Equivalent to: intersectBed -a exons -b rmsk -s
    # Uses:           IntervalFile.all_hits(same_strand=True)
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        # use "same_strand" to enforce, well, same strand.
        for rmsk_hit in rmsk.all_hits(exon, same_strand=True):
            print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, exon.strand,
                                             rmsk_hit.chrom, rmsk_hit.start, rmsk_hit.end, rmsk_hit.strand])



    ##########################################################
    # ex7. Report overlap b/w exons and rmsk where the rmsk 
    #    feature covers at least 50% of the exon.
    #
    # Equivalent to: intersectBed -a exons -b rmsk -f 0.50
    # Uses:           IntervalFile.all_hits()
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        # use "ovlp_pct" to enforce the faction of overlap w.r.t to exon
        for rmsk_hit in rmsk.all_hits(exon, ovlp_pct=0.50):
            print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, exon.strand,
                                             rmsk_hit.chrom, rmsk_hit.start, rmsk_hit.end, rmsk_hit.strand])



    ##########################################################
    # ex8. Report overlap b/w exons and rmsk where the rmsk 
    #    feature covers at least 50% of the exon.
    #
    # Equivalent to: intersectBed -a exons -b rmsk -s -f 0.50
    # Uses:           IntervalFile.all_hits()
    ##########################################################
    exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    for exon in exons:
        # use "same_strand" to enforce, well, same strand.
        for rmsk_hit in rmsk.all_hits(exon, same_strand=True, ovlp_pct=0.50):
            print "\t".join(str(f) for f in [exon.chrom, exon.start, exon.end, exon.strand,
                                             rmsk_hit.chrom, rmsk_hit.start, rmsk_hit.end, rmsk_hit.strand])



if __name__ == "__main__":
    main()
