from libc.stdio cimport FILE
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string as cppstring
from libcpp.memory cimport unique_ptr
from libsequence.polysitevector cimport polymorphicSite,polySiteVector,psite_vec_const_itr,psite_vec_itr

cdef extern from "Sequence/PolyTable.hpp" namespace "Sequence" nogil:
    cdef cppclass CppPolyTable "Sequence::PolyTable":
        const string & operator[](const size_t &) const
        vector[string].iterator begin()
        vector[string].iterator end()
        vector[double].iterator pbegin()
        vector[double].iterator pend()
        vector[pair[double,string]].const_iterator sbegin() const
        vector[pair[double,string]].const_iterator send() const
        vector[string] GetData() const
        vector[double] GetPositions() const

        bint assign(const vector[pair[double,cppstring]].iterator &,const vector[pair[double,cppstring]].iterator & )
        bint assign(vector[double], vector[string])
        bint empty() const
        double position(unsigned &) const
        unsigned numsites() const
        unsigned size() const
        
                
cdef extern from "Sequence/SimData.hpp" namespace "Sequence" nogil:
    cdef cppclass CppSimData "Sequence::SimData"(CppPolyTable):
        CppSimData()
        CppSimData(const CppSimData &)
        CppSimData(const psite_vec_itr,const psite_vec_itr)
        int fromfile(FILE * f)

cdef extern from "Sequence/PolySites.hpp" namespace "Sequence" nogil:
    cdef cppclass CppPolySites "Sequence::PolySites"(CppPolyTable):
        CppPolySites()
        CppPolySites(const psite_vec_itr,const psite_vec_itr)

cdef extern from "Sequence/PolyTableFunctions.hpp" namespace "Sequence" nogil:
    bint polyTableValid(const CppPolyTable * t)
    T cpp_removeGaps"Sequence::removeGaps"[T](const T &, const bint skipAnc, const unsigned anc, const char gapchar)
    T removeInvariantPos[T](const T & t, const bint skipAnc, const unsigned anc,  const char gapchar)
    T cpp_removeAmbiguous"Sequence::removeAmbiguous"[T](const T & t, const bint skipAnc, const unsigned anc, const char gapchar)
    T cpp_removeMissing"Sequence::removeMissing"[T](const T & t, const bint skipAnc, const unsigned anc, const char gapchar)
    T cpp_removeMultiHits"Sequence::removeMultiHits"[T](const T & t, const bint skipAnc, const unsigned anc, const char gapchar)
    T polyTableToBinary[T](const T & t, const unsigned ref , const char gapchar)
    T polyTableFreqFilter[T](const T & t, const unsigned mincount,const bint skipAnc, const unsigned anc, const char gapchar)
    
cdef class PolyTable:
    cdef unique_ptr[CppPolyTable] thisptr
    cpdef size(self)
    cpdef numsites(self)
    cpdef data(self)
    cpdef pos(self)
    cpdef empty(self)
    cpdef assign(self,const vector[polymorphicSite] & d)
    cpdef assign_sep(self,const vector[double] & pos,const vector[string] & data)

cdef class SimData(PolyTable):
    pass

cdef class PolySites(PolyTable):
    pass
