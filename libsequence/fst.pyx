from cython.view cimport array as cvarray
from cpython cimport array
cdef class fst:
    """
    "Factory" object for :math:`F_{st}` calculations.
    """  
    def __cinit__(self, polyTable p,const unsigned[:] config, const double[:] weights = None, bint haveOutgroup = False, unsigned outgroup = 0):
        if weights is None:
            self.thisptr = new FST(p.thisptr,len(config),&config[0],NULL,haveOutgroup,outgroup)
        else:
            if len(config) != len(weights):
                raise RuntimeError("len(config) must equal len(weights)")
            self.thisptr = new FST(p.thisptr,len(config),&config[0],&weights[0],haveOutgroup,outgroup)
    def hsm(self):
        """
        Hudson, Slatkin, Maddison
        """
        return self.thisptr.HSM()
    def slatkin(self):
        """
        Slatkin
        """
        return self.thisptr.Slatkin()
    def hbk(self):
        """
        Hudson, Boos, and Kaplan
        """
        return self.thisptr.HBK()
    def piB(self):
        """
        :math:`\\pi` between populations
        """
        return self.thisptr.piB()
    def piT(self):
        """
        Total :math:`\\pi`
        """
        return self.thisptr.piT()
    def piS(self):
        """
        :math:`\\pi_S`
        """
        return self.thisptr.piS()
    def piD(self):
        """
        :math:`\\pi_D`
        """
        return self.thisptr.piD()
    def shared(self,unsigned i,unsigned j):
        """
        Returns positions of mutations shared between demes i and j
        """
        return self.thisptr.shared(i,j)
    def priv(self,unsigned i,unsigned j):
        """
        Returns set of mutations private to i and private to j, when
        only demes i and j are compared
        """
        return self.thisptr.Private(i,j)
    def fixed(self,unsigned i,unsigned j):
        """
        Returns set of fixed differences between i and j (only considering that comparison)
        """
        return self.thisptr.fixed(i,j)
