from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.set cimport set

from libsequence.polytable cimport polyTable,PolyTable

cdef extern from "Sequence/FST.hpp" namespace "Sequence":
    cdef cppclass FST:
        FST( const PolyTable *, unsigned npop, const unsigned * config,
             const double * weights, bint haveOutgroup, unsigned outgroup ) except +
        double HSM() const
        double Slatkin() const
        double HBK() const
        double piB() const
        double piT() const
        double piS() const
        double piD() const
        set[double] shared(unsigned pop1, unsigned pop2) const
        set[double] fixed(unsigned pop1, unsigned pop2) const
        pair[set[double],set[double]] Private(unsigned pop1, unsigned pop2) const

cdef class fst:
    cdef FST * thisptr
