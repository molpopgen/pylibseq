from libcpp.vector cimport vector
from libsequence.polysitevector cimport polySiteVector,psite_vec_itr,psite_vec_const_itr
from libsequence.polytable cimport CppSimData,PolyTable

cdef class SimDataVec:
    cdef vector[CppSimData] vec
