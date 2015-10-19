# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from cython.operator import dereference as deref,postincrement as inc

from polytable cimport PolyTable,SimData

cdef class polyTable:
    cdef PolyTable * thisptr
    def __cinit__(self):
        thisptr = NULL
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    def size(self):
        return self.thisptr.second.size()
    def __getitem__(self, size_t i):
        return self.thisptr.second[i]
    
cdef class simData(polyTable):
    def __cinit__(self,const vector[polymorphicSite] & data):
        self.thisptr = new SimData(data.const_begin(),data.const_end())
    def __dealloc__(self):
        pass

