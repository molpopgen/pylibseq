from libcpp.vector cimport vector
from libsequence.polysitevector cimport polySiteVector,psite_vec_itr,psite_vec_const_itr
from libsequence.polytable cimport SimData,PolyTable

cdef class simDataVec:
    cdef vector[SimData] vec
