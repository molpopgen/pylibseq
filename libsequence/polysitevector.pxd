from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Sequence/polySiteVector.hpp" namespace "Sequence":
    ctypedef pair[double,string] polymorphicSite
    ctypedef vector[polymorphicSite] polySiteVector
    
ctypedef vector[polymorphicSite].const_iterator psite_vec_const_itr
ctypedef vector[polymorphicSite].iterator psite_vec_itr

        
