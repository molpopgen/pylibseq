# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from cython.operator import dereference as deref,postincrement as inc

from polytable cimport PolyTable,SimData

cdef class polyTable:
    def __cinit__(self):
        thisptr = NULL
    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
    cpdef size(self):
        assert (self.thisptr != NULL)
        return self.thisptr.second.size()
    def __getitem__(self, size_t i):
        assert (self.thisptr != NULL)
        return self.thisptr.second[i]
    cpdef GetData(self):
        assert (self.thisptr != NULL)
        return self.thisptr.GetData()
    cpdef GetPositions(self):
        assert (self.thisptr != NULL)
        return self.thisptr.GetPositions()
    cpdef empty(self):
        assert (self.thisptr != NULL)
        return self.thisptr.empty()
    cpdef assign(self,const vector[polymorphicSite] & d):
        assert (self.thisptr != NULL)
        cdef bint rv = self.thisptr.assign(d.const_begin(),d.const_end())
        if rv == False:
            raise RuntimeError("assign failed")
    cpdef assign_sep(self,const vector[double] & pos,const vector[string] & data):
        assert (self.thisptr != NULL)
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
