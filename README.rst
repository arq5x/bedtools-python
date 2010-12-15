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
  2. GNU compiler

----------------
2. Installation
----------------

::

	python setup.py install
	./test.sh
	# Try it out
	./example-intersect.py
	
---------------
3. Examples
---------------


3a. Intersect one file against another (akin to intersectBed)
---------------------------------------------------------

::

	>>> snps = IntervalFile("testData/snps.hg18.chr21.bed")
	>>> rmsk = IntervalFile("testData/rmsk.hg18.chr21.bed")
	>>>
	>>> # find snps that overlap with repeat annotations
	>>> for a in snps:
	>>>     for hit in rmsk.search(a.chrom, a.start, a.end):
	>>>        print a.chrom, a.start, a.end, a.name,
	>>>        print hit.chrom, hit.start, hit.end, hit.name
	chr21 9719803 9719804 rs55981545 chr21 9719768 9721892 ALR/Alpha
	chr21 9719863 9719864 rs73327798 chr21 9719768 9721892 ALR/Alpha
	chr21 9719950 9719951 rs73327799 chr21 9719768 9721892 ALR/Alpha
	chr21 9719972 9719973 rs71245703 chr21 9719768 9721892 ALR/Alpha
	chr21 9719980 9719981 rs28971396 chr21 9719768 9721892 ALR/Alpha
	...