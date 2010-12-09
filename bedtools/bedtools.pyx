"""
    bedtools.pyx: A Cython wrapper for the BEDTools BedFile class
    
    Author: Aaron Quinlan
    Affl:   Center for Public Health Genomics
            University of Virginia
    Email:  aaronquinlan at gmail dot com
"""
from cython.operator cimport dereference as deref, preincrement as inc #dereference and increment operators
from cpython cimport bool 

"""
    Expose STL strings to Cython.
"""
cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()

from libcpp.vector cimport vector


            
"""
    ***********************************************
    STATUS:  Doesn't yet work.  'Splainin' to do.
    
    GOAL:
        Expose the most basic interface.  Namely,
        load a file into a UCSC-bin tree and search
        for overlaps in the tree.
        
    TO DO:
        1. Need to expose the BED struct.
        2. Need to expose STL vectors.
        3. Need to expose the BedLineStatus enum.
    ***********************************************
"""


     
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
        BedFile(string bedFile)
        Open()
        Close()
        BedLineStatus GetNextBed(BED &bed, int &lineNum)
        loadBedFileIntoMap()

        vector[BED] FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand, bool forceStrand)
        

cdef class Bed:
    cdef BED _bed

cdef Bed create_bed(BED b):
    cdef Bed pyb = Bed.__new__(Bed)
    pyb._bed = b
    return pyb

cdef list vec2list(vector[BED] bv):
    cdef list l = []
    cdef size_t size, i
    cdef BED b
    size = bv.size()
    for i in range(size):
        l.append(create_bed(bv.at(i)))
    return l



        
cdef class IntervalFile:
    """
        Now, the important business: exposing the Interval class to Python code.

        The future.  Something like:
        
        from bedtools-python import IntervalFile

        >>> rmskFile = IntervalFile(rmsk.bed)
        >>> rmskFile.load()
        
        >>> chrom  = "chr1"
        >>> start  = 10
        >>> end    = 20
        >>> strand = "+"
        >>> forceSameStand = False
        >>> hits = rmskFile.findOverlaps(chrom, start, end, strand, forceSameStand)
        
    """
    cdef BedFile *intervalFile_ptr

    def __init__(self, intervalFile):
        self.intervalFile_ptr = new BedFile(string(intervalFile))

    def __dealloc__(self):
        del self.intervalFile_ptr
        
    def loadIntoMap(self):
        self.intervalFile_ptr.loadBedFileIntoMap()
        
    def findOverlaps(self, chrom, start, end, strand, forceStrand):
        """
        ISSUES:
            1. Python will expect hits to be returned as a list, not in a pass-by-reference manner.  FIX.
        """
        return vec2list(self.intervalFile_ptr.FindOverlapsPerBin(string(chrom), int(start), int(end), string(strand), bool(forceStrand)))



