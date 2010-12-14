from cpython cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char *c_str()



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
    cdef cppclass BED:
        string chrom
        CHRPOS start
        CHRPOS end
        string name
        string score
        string strand
        vector[string] otherFields
        vector[BED] overlaps
        BED()
        BED(string chrom, CHRPOS start, CHRPOS end, string name,
             string score, string strand, vector[string] otherFields)


    cdef cppclass BedFile:
        BedFile(string)
        void Open()
        void Close()
        BedLineStatus GetNextBed(BED &bed, int &lineNum)
        void loadBedFileIntoMap()

        # this version doesn't care if the strands match.
        vector[BED] FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, float overlapFraction)
        # if strand is passed, require that the strands match,
        vector[BED] FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand, float overlapFraction)
