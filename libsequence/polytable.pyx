# distutils: language = c++
from libc.stdio cimport * 
from libcpp cimport nullptr
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.cast cimport dynamic_cast
from cython.operator import dereference as deref,postincrement as inc

from polytable cimport CppPolyTable,CppSimData
from polytable cimport polyTableFreqFilter as ptFF

cdef class PolyTable:
    """
    libsequence's Sequence::PolyTable

    .. note:: It is an error to directly use this class in Python.  An assertion will be triggered.  This is simply the base API for other types.
    """
    def __cinit__(self):
        if self.__class__ == PolyTable:
            raise RuntimeError("PolyTable cannot be used directly.  Use derived types instead")
        self.thisptr.reset(nullptr)
    def __dealloc__(self):
        self.thisptr.reset(nullptr)
    cpdef size(self):
        """
        Get the sample size of the polymorphism table
        """
        assert (self.thisptr != NULL)
        return self.thisptr.get().size()
    def __getitem__(self, size_t i):
        #assert (self.thisptr != NULL)
        return deref(self.thisptr)[i]
    cpdef numsites(self):
        return self.thisptr.get().numsites()
    cpdef data(self):
        """
        Get the genotype data

        :rtype: list
        """
        #assert (self.thisptr != NULL)
        return self.thisptr.get().GetData()
    cpdef pos(self):
        """
        Get the mutation positions

        :rtype: list
        """
        #assert (self.thisptr != NULL)
        return self.thisptr.get().GetPositions()
    cpdef empty(self):
        """
        Return True if sample size is 0
        """
        #assert (self.thisptr != NULL)
        return self.thisptr.get().empty()
    cpdef assign(self,const vector[polymorphicSite] & d):
        """
        Fill data from a list of tuples

        Example:

        >>> import libsequence.polytable as pypt
        >>> x = pypt.SimData()
        >>> x.assign([ (0.1,"01"),(0.2,"10") ])
        """
        #assert (self.thisptr != NULL)
        cdef bint rv = self.thisptr.get().assign(d.const_begin(),d.const_end())
        if rv == False:
            raise RuntimeError("assign failed")
    cpdef assign_sep(self,const vector[double] & pos,const vector[string] & data):
        """
        Fill object from two lists

        Example: 

        >>> import libsequence.polytable as pypt
        >>> x = pypt.SimData()
        >>> pos = [0.1,0.2,0.3,0.4]
        >>> data = ["0101","1011"]
        >>> x.assign_sep(pos,data)
        """
        assert (self.thisptr != NULL)
        #cdef bint rv =self.thisptr.assign[double,string](pos.data(),pos.size(),data.data(),data.size())
        cdef bint rv =self.thisptr.get().assign(pos,data)
        if rv == False:
            raise RuntimeError("assign_sep failed")
    def tolist(self):
        """
        Return data as list of tuples.
        """
        cdef vector[pair[double,string]] rv
        rv.assign(self.thisptr.get().sbegin(),self.thisptr.get().send())
        return rv
    
cdef class SimData(PolyTable):
    """
    A polymorphism table for binary data.  0/1 = ancestral/derived.

    .. note:: See :class:`libsequence.polytable.PolyTable`
    """
    def __cinit__(self):
        self.thisptr = <unique_ptr[CppPolyTable]>unique_ptr[CppSimData](new CppSimData())
    def __init__(self,object x = None, object y = None):
        """
        Constructor

        :param x: A list
        :param y: A list

        If x is None and y is None, and empty object is created.
        
        If y is None, x must be a list of tuples (position, "string")

        if neither x nor y are None, x is a list of positions, and y a list of haplotypes.

        Example:

        >>> import libsequence.polytable as polyt
        >>> x = polyt.SimData()
        >>> x = polyt.SimData( [ (0.1,"0101"),(0.2,"1010") ] )
        >>> x = polyt.SimData( [0.2,0.2],["01","10","01","10"] )
        """
        if x is None:
            if y is None:
                pass
            else:
                raise RuntimeError("y must be None if x is None")
        else:
            if y is None:
                self.assign(x)
            else:
                self.assign_sep(x,y)

cdef class PolySites(PolyTable):
    """
    A polymorphism table for Sequence data.  A/G/C/T/N/0/1 is the alphabet for analysis.

    .. note:: See :class:`libsequence.Polytable.PolyTable`
    """
    def __cinit__(self):
        self.thisptr = <unique_ptr[CppPolyTable]>unique_ptr[CppPolySites](new CppPolySites())
    def __init__(self,list x = None, list y = None):
        """
        Constructor

        :param x: A list
        :param y: A list

        If x is None and y is None, and empty object is created.
        
        If y is None, x must be a list of tuples (position, "string")

        if neither x nor y are None, x is a list of positions, and y a list of haplotypes.

        Example:

        >>> import libsequence.polytable as polyt
        >>> x = polyt.PolySites()
        >>> x = polyt.PolySites( [ (0.1,"ATTA"),(0.2,"GGGC") ] )
        >>> x = polyt.PolySites( [0.2,0.2],["AG","TG","TG","AC"] )

        """
        if x is None:
            if y is None:
                pass
            else:
                raise RuntimeError("y must be None if x is None")
        else:
            if y is None:
                self.assign(x)
            else:
                self.assign_sep(x,y)
                
def removeGaps(PolyTable p, gapchar = '-'):
    """
    Remove all sites (columns) with gaps.

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    :param gapchar: the character representing an alignment gap
    
    Example:

    >>> import libsequence.polytable as pypt
    >>> x = pypt.PolySites()
    >>> x.assign_sep([0.1,0.2,0.3],["ATC","CGA","AT-"])
    >>> pypt.removeGaps(x)
    >>> x.pos()
    [0.1, 0.2]
    """
    cdef bytes gc = <bytes>gapchar
    cdef CppPolySites temp
    cdef CppSimData temp2
    if isinstance(p,PolySites):
        temp = cpp_removeGaps[CppPolySites](deref(dynamic_cast['CppPolySites*'](p.thisptr.get())),False,0,gc[0])
        p.assign_sep(temp.GetPositions(),temp.GetData())
    else:
        temp2 = cpp_removeGaps[CppSimData](deref(dynamic_cast['CppSimData*'](p.thisptr.get())),False,0,gc[0])
        p.assign_sep(temp2.GetPositions(),temp2.GetData())

def isValid(PolyTable p):
    """
    Return True if p only contains the characters A,C,G,T,N,-,0,1.

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    
    .. note:: This is not case-sensitive

    Example:

    >>> import libsequence.polytable as pypt
    >>> x = pypt.PolySites()
    >>> x.assign_sep([0.1,0.2,0.3],["ATC","CGA","AT-"])
    >>> pypt.isValid(x)
    True

    Now, Z is not a valid character:
    
    >>> x.assign_sep([0.1,0.2,0.3],["ATC","CGA","ATZ"])
    >>> pypt.isValid(x)
    False
    """
    return polyTableValid(p.thisptr.get())

def removeMono(PolyTable p, bint skipOutgroup = False, unsigned outgroup = 0, gapchar = '-'):
    """
    Remove invariant sites from p

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    :param skipOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    cdef bytes gc = <bytes>gapchar
    cdef CppPolySites temp
    cdef CppSimData temp2
    if isinstance(p,PolySites):
        temp = removeInvariantPos[CppPolySites](deref(dynamic_cast['CppPolySites*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp.GetPositions(),temp.GetData())
    else:
        temp2 = removeInvariantPos[CppSimData](deref(dynamic_cast['CppSimData*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp2.GetPositions(),temp2.GetData())
        
def freqFilter(PolyTable p,unsigned mincount,bint skipOutgroup = False, unsigned outgroup = 0,gapchar='-'):
    """
    Remove all sites with mutation count <= mincount

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    :param mincount: Sites with mutation counts <= this value will be removed.  For example, use a value of 2 to remove singletons
    :param haveOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    
    .. note:: If haveOutgroup == True, this is a filter on derived mutation counts, otherwise it is a filter on minor allele counts.  If p is of type :class:`libsequence.PolyTable.SimData`, this is a filter on derived mutation counts.
    """
    cdef bytes gc = <bytes>gapchar
    cdef CppPolySites temp
    cdef CppSimData temp2
    if isinstance(p,PolySites):
        temp = ptFF[CppPolySites](deref(dynamic_cast['CppPolySites*'](p.thisptr.get())),mincount,skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp.GetPositions(),temp.GetData())
    else:
        temp2 = ptFF[CppSimData](deref(dynamic_cast['CppSimData*'](p.thisptr.get())),mincount,skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp2.GetPositions(),temp2.GetData())

def removeMissing(PolyTable p,bint skipOutgroup = False, unsigned outgroup = 0,gapchar='-'):
    """
    Remove all sites missing data (the 'N' or 'n' character)

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    :param mincount: Sites with mutation counts <= this value will be removed.  For example, use a value of 2 to remove singletons
    :param skipOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    cdef bytes gc = <bytes>gapchar
    cdef CppPolySites temp
    cdef CppSimData temp2
    if isinstance(p,PolySites):
        temp = cpp_removeMissing[CppPolySites](deref(dynamic_cast['CppPolySites*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp.GetPositions(),temp.GetData())
    else:
        temp2 = cpp_removeMissing[CppSimData](deref(dynamic_cast['CppSimData*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp2.GetPositions(),temp2.GetData())
        
def removeMultiHits(PolyTable p,bint skipOutgroup = False, unsigned outgroup = 0,gapchar='-'):
    """
    Remove all sites with more than two character states.

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    :param skipOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    cdef bytes gc = <bytes>gapchar
    cdef CppPolySites temp
    cdef CppSimData temp2
    if isinstance(p,PolySites):
        temp = cpp_removeMultiHits[CppPolySites](deref(dynamic_cast['CppPolySites*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp.GetPositions(),temp.GetData())
    else:
        temp2 = cpp_removeMultiHits[CppSimData](deref(dynamic_cast['CppSimData*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp2.GetPositions(),temp2.GetData())
    
def removeAmbiguous(PolyTable p,bint skipOutgroup = False, unsigned outgroup = 0,gapchar='-'):
    """
    Remove all sites with characters not in the set A,G,C,T,N,0,1,-.

    :param p: An object derived from :class:`libsequence.polytable.PolyTable`
    :param skipOutgroup: if True, the sequence at position 'outgroup' will not be included in determining if a site is monomorphic
    :param outgroup: The index of the outgroup sequence in p
    """
    cdef bytes gc = <bytes>gapchar
    cdef CppPolySites temp
    cdef CppSimData temp2
    if isinstance(p,PolySites):
        temp = cpp_removeAmbiguous[CppPolySites](deref(dynamic_cast['CppPolySites*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp.GetPositions(),temp.GetData())
    else:
        temp2 = cpp_removeAmbiguous[CppSimData](deref(dynamic_cast['CppSimData*'](p.thisptr.get())),skipOutgroup,outgroup,gc[0])
        p.assign_sep(temp2.GetPositions(),temp2.GetData())

def readSimData(SimData d):
    """
    Read a :class:`libsequence.polytable.SimData` from the standard
    input stream.

    Intended use would be writing a program to read data in from a pipe.

    .. code-block:: python

        import libsequence.polytable as pt
        d=pt.SimData()
        while pt.readSimData(d) is True:
            print d.numsites()
    
    :return: True if object could be read, False if stream is empty
    """
    if feof(stdin):
        return False
    (dynamic_cast['CppSimData*'](d.thisptr.get())).fromfile(stdin)
    return True
