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
    examples/ex0-echo.py
    examples/ex1-intersect.py
    examples/ex2-pysam-and-bedtools.py
    examples/ex3-fasta-from-bed.py
    examples/ex4-window.py

----------------
3. The API
----------------
Still need to document this.

---------------
4. Examples
---------------

Below is a brief example to whet your appetite.  For more details on how to use the API, check out the examples/ directory.
::

  >>> from bedtools import IntervalFile
  >>>
  >>>
  >>> exons = IntervalFile("bedtools/tests/data/exons.hg18.chr21.bed")
  >>> rmsk  = IntervalFile("bedtools/tests/data/rmsk.hg18.chr21.bed")
  >>>
  >>> # find exons that overlap with repeat annotations
  >>> ex = exons.next()
  >>> for hit in rmsk.search(ex):
  ...      print ex.chrom,  ex.start,  ex.end,\
  ...            hit.chrom, hit.start, hit.end,
  chr21 9928613 9928911 chr21 9928614 9928678
