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

---------------
2. Purpose
---------------

---------------
3. Limitations
---------------

---------------
4. Examples
---------------

Load A Bed File
---------------

::

    >>> from bedtools import IntervalFile
    >>> bed = IntervalFile('bedtools/tests/data/rmsk.hg18.chr21.bed')
    >>> bed.loadIntoMap()
    >>> bed.findOverlaps("chr21", 9719768, 9724768)
    [Bed(chr21:9719768..9721892), Bed(chr21:9721905..9725582)]
