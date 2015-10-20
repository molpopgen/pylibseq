from cython.view cimport array as cvarray
from cpython cimport array
cdef class fst:
    def __cinit__(self, polyTable p,const unsigned[:] config, const double[:] weights = None, bint haveOutgroup = False, unsigned outgroup = 0):
        if weights is None:
            self.thisptr = new FST(p.thisptr,len(config),&config[0],NULL,haveOutgroup,outgroup)
        else:
            if len(config) != len(weights):
                raise RuntimeError("len(config) must equal len(weights)")
            self.thisptr = new FST(p.thisptr,len(config),&config[0],&weights[0],haveOutgroup,outgroup)
    def hsm(self):
        return self.thisptr.HSM()
    def slatkin(self):
        return self.thisptr.Slatkin()
    def hbk(self):
        return self.thisptr.HBK()
    def piB(self):
        return self.thisptr.piB()
    def piT(self):
        return self.thisptr.piT()
    def piS(self):
        return self.thisptr.piS()
    def piD(self):
        return self.thisptr.piD()
    def shared(self,unsigned i,unsigned j):
        return self.thisptr.shared(i,j)
    def priv(self,unsigned i,unsigned j):
        return self.thisptr.Private(i,j)
    def fixed(self,unsigned i,unsigned j):
        return self.thisptr.fixed(i,j)
