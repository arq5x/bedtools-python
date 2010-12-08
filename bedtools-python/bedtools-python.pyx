"""
    BedFile.pyx: A Cython wrapper for the BEDTools BedFile class
    
    Author: Aaron Quinlan
    Affl:   Center for Public Health Genomics
            University of Virginia
    Email:  aaronquinlan at gmail dot com
"""



"""
    Expose STL strings to Cython.
"""
cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()


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
    cdef cppclass BedFile:
        BedFile(string &)
        Open()
        Close()
        BedLineStatus GetNextBed(BED &bed, int &lineNum);
        loadBedFileIntoMap();

        FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand, vector<BED> &hits, bool forceStrand);
        
        
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

    def __init__(self, chrom, start, end):
        self.intervalFile_ptr = new Interval(string(intervalFile))

    def __dealloc__(self):
        del self.intervalFile_ptr
        
    def load(self):
        intervalFile_ptr.loadBedFileIntoMap()
        
    def findOverlaps(self, chrom, start, end, strand, hits, forceStrand):
        """
        ISSUES:
            1. Python will expect hits to be returned as a list, not in a pass-by-reference manner.  FIX.
        """
        intervalFile_ptr.FindOverlapsPerBin()
    


