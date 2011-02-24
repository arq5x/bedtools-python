:Project: bedtools-python
:Version: 0.1.0
:Authors: - Aaron Quinlan, University of Virginia
          - Brent Pedersen, University of California, Berkeley
:Contact: arq5x@virginia.edu

===============
bedtools-python
===============

---------------
1. Requirements
---------------
  1. Cython
    - sudo easy_install cython
  2. GNU compiler
  3. zlib

----------------
2. Installation
----------------

::

    python setup.py install
    ./test.sh
    # Try it out
    ./example-intersect.py

----------------
3. The API
----------------
Still need to document this.

---------------
4. Examples
---------------
4a. Intersect one file against another (akin to intersectBed)
--------------------------------------------------------------

::

    >>> from bedtools import IntervalFile, Overlap
    >>> 
    >>> def main():
    >>> 
    >>> exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
    >>> rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
    >>> 
    >>> # find exons that overlap with repeat annotations
    >>> for ex in exons:
    >>>     # retrieve repeats that overlap this exon
    >>>     for hit in rmsk.search(ex.chrom, ex.start, ex.end):
    >>>         # report the _full_ features that overlap
    >>>         print ex.chrom,  ex.start,  ex.end,\
    >>>               hit.chrom, hit.start, hit.end,
    >>>         # extract the coordinates of the overlap, 
    >>>         # as well as the amount of overlap in b.p.
    >>>         o = Overlap(ex.start, ex.end, hit.start, hit.end)
    >>>         print "overlap_start=" + str(o.overlap_start),\
    >>>               " overlap_end="  + str(o.overlap_end),\
    >>>               " overlap_amt="  + str(o.overlap_amt)

    chr21 10119416 10119517 chr21 10119432 10119526 overlap_start=10119432  overlap_end=10119517  overlap_amt=85
    chr21 10120573 10120796 chr21 10119986 10121803 overlap_start=10120573  overlap_end=10120796  overlap_amt=223
    chr21 10119416 10119517 chr21 10119432 10119526 overlap_start=10119432  overlap_end=10119517  overlap_amt=85
    chr21 10120593 10120808 chr21 10119986 10121803 overlap_start=10120593  overlap_end=10120808  overlap_amt=215
    chr21 10081596 10081687 chr21 10081575 10081768 overlap_start=10081596  overlap_end=10081687  overlap_amt=91
    
    
    
    
    
    ...