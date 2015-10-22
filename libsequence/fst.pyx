from cython.view cimport array as cvarray
from cpython cimport array
cdef class fst:
    """
    "Factory" object for :math:`F_{st}` calculations.
    """  
    def __cinit__(self, polyTable p,list config, list weights = None, bint haveOutgroup = False, unsigned outgroup = 0):
        c = array.array('I',config)
        cdef unsigned[:] cv = c
        cdef double[:] wv
        if weights is None:
            self.thisptr = new FST(p.thisptr,len(config),&cv[0],NULL,haveOutgroup,outgroup)
        else:
            w = array.array('d',weights)
            wv = w
            if len(config) != len(weights):
                raise RuntimeError("len(config) must equal len(weights)")
            self.thisptr = new FST(p.thisptr,len(config),&cv[0],&wv[0],haveOutgroup,outgroup)
    def __dealloc__(self):
        del self.thisptr
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

        :rtype: list
        
        """
        return sorted(list(self.thisptr.shared(i,j)))
    def priv(self,unsigned i,unsigned j):
        """
        Returns set of mutations private to i and private to j, when
        only demes i and j are compared

        :return: dict. i and j are keys, and list of private positions are values
        
        :rtype: dict
        """
        p = self.thisptr.Private(i,j)
        return {i:sorted(list(p.first)),j:sorted(list(p.second))}
    def fixed(self,unsigned i,unsigned j):
        """
        Returns set of fixed differences between i and j (only considering that comparison)

        :rtype: list
        """
        f = self.thisptr.fixed(i,j)
        return sorted(list(f))
