"""
    bedtools.pyx: A Cython wrapper for the BEDTools BedFile class
    
    Author: Aaron Quinlan
    Affl:   Center for Public Health Genomics
            University of Virginia
    Email:  aaronquinlan at gmail dot com
"""
include "cbedtools.pxi"
from cython.operator cimport dereference as deref

cdef class Bed:
    cdef BED *_bed

    @property
    def start(self):
        return self._bed.start

    @property
    def chrom(self):
        return self._bed.chrom.c_str()

    @property
    def end(self):
        return self._bed.end

    @property
    def name(self):
        return self._bed.name.c_str()

    def __repr__(self):
        return "Bed(%s:%i..%i)" % (self._bed.chrom.c_str(), self._bed.start, self._bed.end)

    def __dealloc__(self):
        del self._bed

cdef Bed create_bed(BED b):
    cdef Bed pyb = Bed.__new__(Bed)
    pyb._bed = new BED(b.chrom, b.start, b.end, b.name, b.score, b.strand, b.otherFields)
    return pyb

cdef list vec2list(vector[BED] bv):
    cdef list l = []
    cdef size_t size = bv.size(), i
    cdef BED b
    for i in range(size):
        b = bv.at(i)
        l.append(create_bed(b))
    return l

cdef class IntervalFile:
    cdef BedFile *intervalFile_ptr

    def __init__(self, intervalFile):
        self.intervalFile_ptr = new BedFile(string(intervalFile))

    def loadIntoMap(self):
        self.intervalFile_ptr.loadBedFileIntoMap()

    def __dealloc__(self):
        del self.intervalFile_ptr

    def findOverlaps(self, chrom, int start, int end, strand, bool forceStrand):
        cdef vector[BED] vec_b = self.intervalFile_ptr.FindOverlapsPerBin(string(chrom), start, end, string(strand), bool(forceStrand))
        try:
            return vec2list(vec_b)
        finally:
            pass
