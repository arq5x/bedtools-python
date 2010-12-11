"""
    bedtools.pyx: A Cython wrapper for the BEDTools BedFile class
    
    Author: Aaron Quinlan
    Affl:   Center for Public Health Genomics
            University of Virginia
    Email:  aaronquinlan at gmail dot com
"""
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
from cpython cimport bool 
from libcpp.vector cimport vector

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()



"""
Create Cython definitions for the Interval API defined in Interval.h
"""
cdef extern from "bedFile.h":
    cdef enum BedLineStatus:
        BED_INVALID = -1
        BED_HEADER  = 0
        BED_BLANK   = 1
        BED_VALID   = 2

    ctypedef unsigned int CHRPOS
    cdef struct BED:
        string chrom
        CHRPOS start
        CHRPOS end
        string name
        string score
        string strand
        vector[string] otherFields
        vector[BED] overlaps

    cdef cppclass BedFile:
        BedFile(string)
        void Open()
        void Close()
        BedLineStatus GetNextBed(BED &bed, int &lineNum)
        void loadBedFileIntoMap()

        vector[BED] FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand, bool forceStrand)


cdef class Bed:
    cdef BED *_bed

cdef Bed create_bed(BED b):
    cdef Bed pyb = Bed.__new__(Bed)
    pyb._bed = &b
    return pyb

cdef list vec2list(vector[BED] bv):
    cdef list l = []
    cdef size_t size = bv.size(), i
    cdef BED b
    for i in range(size):
        l.append(create_bed(bv.at(i)))
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
        return vec2list(self.intervalFile_ptr.FindOverlapsPerBin(string(chrom), start, end, string(strand), bool(forceStrand)))
