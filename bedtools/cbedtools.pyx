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
    def chrom(self):
        return self._bed.chrom.c_str()

    @property
    def start(self):
        return self._bed.start
            
    @property
    def end(self):
        return self._bed.end
    
    @property
    def name(self):
        return self._bed.name.c_str()
        
    @property
    def strand(self):
        return self._bed.strand.c_str()
        
    @property
    def other(self):
        return string_vec2list(self._bed.otherFields)
            
    def __repr__(self):
        return "Bed(%s:%i..%i)" % (self._bed.chrom.c_str(), self._bed.start, self._bed.end)

    def __dealloc__(self):
        del self._bed

cdef Bed create_bed(BED b):
    cdef Bed pyb = Bed.__new__(Bed)
    pyb._bed = new BED(b.chrom, b.start, b.end, b.name, b.score, b.strand, b.otherFields)
    return pyb

cdef list string_vec2list(vector[string] sv):
    cdef size_t size = sv.size(), i
    cdef list l = []
    for i in range(size):
        l.append(sv.at(i).c_str())
    return l

cdef list bed_vec2list(vector[BED] bv):
    cdef size_t size = bv.size(), i
    cdef list l = []
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

    def findOverlaps(self, chrom, int start, int end, strand = None, float overlapFraction = 0.0):
        """
        If strand is not passed, hits will be reported without regard to strand.
        If strand is passed, hits will only be reported if on the same strand.
        
        The overlapFraction defaults to 0.0 so that even 1bp of overlap is sufficient.
        If overlapFraction is passed, this fraction of the passed BED feature must be "covered."
        
        Examples:
        1. Find all overlaps regardless of degree on either strand
        findOverlaps("chr1", 10, 20)

        2. Find all overlaps regardless of degree on positive strand
        findOverlaps("chr1", 10, 20, "+")

        3. Find all overlaps covering at least half of this feature on either strand
        findOverlaps("chr1", 10, 20, None, 0.5)
                
        3. Find all overlaps covering at least half of this feature on negative strand
        findOverlaps("chr1", 10, 20, "-", 0.5)
        
        """
        cdef vector[BED] vec_b
        
        if strand is None:
            vec_b = self.intervalFile_ptr.FindOverlapsPerBin(string(chrom), start, end, overlapFraction)
            try:
                return bed_vec2list(vec_b)
            finally:
                pass
        else:
            vec_b = self.intervalFile_ptr.FindOverlapsPerBin(string(chrom), start, end, string(strand), overlapFraction)
            try:
                return bed_vec2list(vec_b)
            finally:
                pass  
