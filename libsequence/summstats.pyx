from libcpp.cast cimport dynamic_cast

cdef class polySNP:
    cdef PolySNP * thisptr
    def __cinit__(self,polyTable p):
       self.thisptr = new PolySNP(p.thisptr,False,0,False)

cdef class polySIM:
    def __cinit__(self,simData d):
        self.thisptr = new PolySIM(dynamic_cast['SimData*'](d.thisptr))
        print self.thisptr.ThetaPi()
    
