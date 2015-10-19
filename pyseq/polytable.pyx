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
    cpdef GetData(self):
        return self.thisptr.GetData()
    cpdef GetPositions(self):
        return self.thisptr.GetPositions()
    cpdef empty(self):
        return self.thisptr.empty()
    def assign(self,const vector[polymorphicSite] & d):
        cdef bint rv = self.thisptr.assign(d.const_begin(),d.const_end())
        if rv == False:
            raise RuntimeError("assign failed")
    def assign_sep(self,const vector[double] & pos,const vector[string] & data):
        cdef bint rv =self.thisptr.assign[double,string](pos.data(),pos.size(),data.data(),data.size())
        if rv == False:
            raise RuntimeError("assign_sep failed")

cdef class simData(polyTable):
    def __cinit__(self):
        self.thisptr = new SimData()
    def __dealloc__(self):
        pass

cdef class polySites(polyTable):
    def __cinit__(self):
        self.thisptr = new PolySites()
    def __dealloc__(self):
        pass
