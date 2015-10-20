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
    cpdef numsites(self):
        return self.thisptr.numsites()
    cpdef data(self):
        """
        Get the genotype data

        :rtype: list
        """
        assert (self.thisptr != NULL)
        return self.thisptr.GetData()
    cpdef pos(self):
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
        >>> x.assign([ (0.1,"01"),(0.2,"10") ])
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
    A polymorphism table for Sequence data.  A/G/C/T/N/0/1 is the alphabet for analysis.

    .. note:: See :class:`libsequence.polytable.polyTable`
    """
    def __cinit__(self):
        self.thisptr = new PolySites()
    def __dealloc__(self):
        pass

def removeGaps(polyTable p, gapchar = '-'):
    """
    Remove all sites (columns) with gaps.

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    :param gapchar: the character representing an alignment gap
    
    Example:

    >>> import libsequence.polytable as pypt
    >>> x = pypt.polySites()
    >>> x.assign_sep([0.1,0.2,0.3],["ATC","CGA","AT-"])
    >>> pypt.removeGaps(x)
    >>> x.pos()
    [0.1, 0.2]
    """
    cdef char * gc = gapchar
    RemoveGaps(p.thisptr,gc[0])

def isValid(polyTable p):
    """
    Return True if p only contains the characters A,C,G,T,N,-,0,1.

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    
    .. note:: This is not case-sensitive

    Example:

    >>> import libsequence.polytable as pypt
    >>> x = pypt.polySites()
    >>> x.assign_sep([0.1,0.2,0.3],["ATC","CGA","AT-"])
    >>> pypt.isValid(x)
    True

    Now, Z is not a valid character:
    
    >>> x.assign_sep([0.1,0.2,0.3],["ATC","CGA","ATZ"])
    >>> pypt.isValid(x)
    False
    """
    return PolyTableValid(p.thisptr)

def removeMono(polyTable p, bint skipOutgroup = False, unsigned outgroup = 0):
    """
    Remove invariant sites from p

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    :param skipOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    RemoveInvariantColumns(p.thisptr,skipOutgroup,outgroup)

def freqFilter(polyTable p,unsigned mincount,bint haveOutgroup = False, unsigned outgroup = 0):
    """
    Remove all sites with mutation count <= mincount

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    :param mincount: Sites with mutation counts <= this value will be removed.  For example, use a value of 2 to remove singletons
    :param haveOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    
    .. note:: If haveOutgroup == True, this is a filter on derived mutation counts, otherwise it is a filter on minor allele counts.  If p is of type :class:`libsequence.polyTable.simData`, this is a filter on derived mutation counts.
    """
    p.thisptr.ApplyFreqFilter(mincount,haveOutgroup,outgroup)

def removeMissing(polyTable p,bint haveOutgroup = False, unsigned outgroup = 0):
    """
    Remove all sites missing data (the 'N' or 'n' character)

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    :param mincount: Sites with mutation counts <= this value will be removed.  For example, use a value of 2 to remove singletons
    :param haveOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    p.thisptr.RemoveMissing(haveOutgroup,outgroup)

def removeMultiHits(polyTable p,bint haveOutgroup = False, unsigned outgroup = 0):
    """
    Remove all sites with more than two character states.

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    :param haveOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    p.thisptr.RemoveMultiHits(haveOutgroup,outgroup)

def removeAmbiguous(polyTable p,bint haveOutgroup = False, unsigned outgroup = 0):
    """
    Remove all sites with characters not in the set A,G,C,T,N,0,1,-.

    :param p: An object derived from :class:`libsequence.polytable.polyTable`
    :param haveOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    p.thisptr.RemoveAmbiguous(haveOutgroup,outgroup)
