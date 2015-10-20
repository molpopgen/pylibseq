from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

from polysitevector cimport polymorphicSite,polySiteVector,psite_vec_const_itr

cdef extern from "Sequence/PolyTable.hpp" namespace "Sequence":
    cdef cppclass PolyTable(pair[vector[double],vector[string]]):
        const string & operator[](const size_t &) const
        vector[string].iterator begin()
        vector[string].iterator end()
        vector[double].iterator pbegin()
        vector[double].iterator pend()
        vector[pair[double,string]].const_iterator sbegin() const
        vector[pair[double,string]].const_iterator send() const
        vector[string] GetData() const
        vector[double] GetPositions() const

        bint assign(const psite_vec_const_itr &,
                    const psite_vec_const_itr & )
        bint assign[NUMERIC,STRING]( const NUMERIC *,
                                     const size_t &,
                                     const STRING *,
                                     const size_t & )
        bint empty() const
        double position(unsigned &) const
        unsigned numsites() const
        unsigned size() const

        #non-const operations
        void ApplyFreqFilter(const unsigned & mincount,
                             const bint & haveOutgroup,
                             const unsigned & outgroup)
        void RemoveMultiHits(const bint & skipOutgroup,
                             const unsigned & outgroup)
        void RemoveMissing(const bint & skipOutgroup,
                           const unsigned & outgroup)
        void RemoveAmbiguous(const bint & skipOutgroup,
                             const unsigned & outgroup)
        
                
cdef extern from "Sequence/SimData.hpp" namespace "Sequence":
    cdef cppclass SimData(PolyTable):
        SimData()
        SimData(const SimData &)

cdef extern from "Sequence/PolySites.hpp" namespace "Sequence":
    cdef cppclass PolySites(PolyTable):
        PolySites()

cdef extern from "Sequence/PolyTableFunctions.hpp" namespace "Sequence":
  void RemoveGaps(PolyTable *t, const char & gapchar)
  void RemoveInvariantColumns(PolyTable *t,
			      const bint & skipOutgroup,
			      const unsigned & outgroup)
  bint PolyTableValid(const PolyTable * t)

cdef class polyTable:
    cdef PolyTable * thisptr
    cpdef size(self)
    cpdef numsites(self)
    cpdef data(self)
    cpdef pos(self)
    cpdef empty(self)
    cpdef assign(self,const vector[polymorphicSite] & d)
    cpdef assign_sep(self,const vector[double] & pos,const vector[string] & data)

cdef class simData(polyTable):
    pass

cdef class polySites(polyTable):
    pass
