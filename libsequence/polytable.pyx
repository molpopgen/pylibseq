# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from cython.operator import dereference as deref,postincrement as inc

from polytable cimport PolyTable,SimData

cdef class polyTable:
    """
    libsequence's Sequence::PolyTable

    .. note:: It is an error to directly use this class in Python.  An assertion will be triggered.  This is simply the base API for other types.
    """
    def __cinit__(self):
        if self.__class__ == polyTable:
            raise RuntimeError("polyTable cannot be used directly.  Use derived types instead")
        thisptr = NULL
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    cpdef size(self):
        """
        Get the sample size of the polymorphism table
        """
        assert (self.thisptr != NULL)
        return self.thisptr.second.size()
    def __getitem__(self, size_t i):
        assert (self.thisptr != NULL)
        return self.thisptr.second[i]
    cpdef GetData(self):
        """
        Get the genotype data

        :rtype: list
        """
        assert (self.thisptr != NULL)
        return self.thisptr.GetData()
    cpdef GetPositions(self):
        """
        Get the mutation positions

        :rtype: list
        """
        assert (self.thisptr != NULL)
        return self.thisptr.GetPositions()
    cpdef empty(self):
        """
        Return True if sample size is 0
        """
        assert (self.thisptr != NULL)
        return self.thisptr.empty()
    cpdef assign(self,const vector[polymorphicSite] & d):
        """
        Fill data from a list of tuples

        Example:

        >>> import libsequence.polytable as pypt
        >>> x = pypt.simData()
        >>>> x.assign([ (0.1,"01"),(0.2,"10") ])
        """
        assert (self.thisptr != NULL)
        cdef bint rv = self.thisptr.assign(d.const_begin(),d.const_end())
        if rv == False:
            raise RuntimeError("assign failed")
    cpdef assign_sep(self,const vector[double] & pos,const vector[string] & data):
        """
        Fill object from two lists

        Example: 

        >>> import libsequence.polytable as pypt
        >>> x = pypt.simData()
        >>> pos = [0.1,0.2,0.3,0.4]
        >>> data = ["0101","1011"]
        >>> x.assign_sep(pos,data)
        """
        assert (self.thisptr != NULL)
        cdef bint rv =self.thisptr.assign[double,string](pos.data(),pos.size(),data.data(),data.size())
        if rv == False:
            raise RuntimeError("assign_sep failed")

cdef class simData(polyTable):
    """
    A polymorphism table for binary data.  0/1 = ancestral/derived.

    .. note:: See :class:`libsequence.polytable.polyTable`
    """
    def __cinit__(self):
        self.thisptr = new SimData()
    def __dealloc__(self):
        pass

cdef class polySites(polyTable):
    """
    A polymorphism table for Sequence data.  0/1 = ancestral/derived.

    .. note:: See :class:`libsequence.polytable.polyTable`
    """
    def __cinit__(self):
        self.thisptr = new PolySites()
    def __dealloc__(self):
        pass

def removeGaps(polyTable p, gapchar = '-'):
    RemoveGaps(p.thisptr,gapchar)

def isValid(polyTable p):
    return PolyTableValid(p.thisptr)

def removeMono(polyTable p, bint skipOutgroup = False, unsigned outgroup = 0):
    RemoveInvariantColumns(p.thisptr,skipOutgroup,outgroup)

def freqFilter(polyTable p,unsigned mincount,bint haveOutgroup = False, unsigned outgroup = 0):
    p.thisptr.ApplyFreqFilter(mincount,haveOutgroup,outgroup)

def removeMissing(polyTable p,bint haveOutgroup = False, unsigned outgroup = 0):
    p.thisptr.RemoveMissing(haveOutgroup,outgroup)

def removeMultiHits(polyTable p,bint haveOutgroup = False, unsigned outgroup = 0):
    p.thisptr.RemoveMultiHits(haveOutgroup,outgroup)

def removeAmbiguous(polyTable p,bint haveOutgroup = False, unsigned outgroup = 0):
    p.thisptr.RemoveAmbiguous(haveOutgroup,outgroup)
