from libcpp.cast cimport dynamic_cast
from cython.view cimport array as cvarray
from cpython cimport array
from cython.operator cimport dereference as deref

cdef extern from "<limits>" namespace "std" nogil:
    cdef cppclass numeric_limits[T]:
        T max()
        T min()
        T epsilon()
        
cdef class polySNP:
    """
    "Factory" class to calculate summary statistics from nucleotide data.

    This is a wrapper around libsequence's Sequence::PolySNP
    """
    def __cinit__(self,polyTable p,bint haveOutgroup = False, unsigned outgroup = 0, bint totMuts = True):
       self.thisptr = new PolySNP(p.thisptr,haveOutgroup,outgroup,totMuts)
    def __dealloc__(self):
        del self.thisptr
    def thetapi(self):
        """
        "Sum of site heterozygosity."  :math:`\\hat\\theta_\\pi = \\sum_i^S\\frac{c}{n}\\frac{n-c-1}{n-1}`,
        where :math:`S` is the number of polymorphisms and :math:`n` is the sample size.
        """
        return self.thisptr.ThetaPi()
    def thetaw (self):
        """
        Watterson's estimator of :math:`\\theta` from :math:`S`
        """
        return self.thisptr.ThetaW()
    def thetah (self):
        """
        Fay and Wu's estimator of :math:`\\theta`.
        """
        return self.thisptr.ThetaH()
    def thetal (self):
        """
        Normalized version of Fay and Wu's estimator of :math:`\\theta`.
        """
        return self.thisptr.ThetaL()

    #calculate various numbers related to polymorphism
    def numpoly(self):
        """
        Number of poylmorphic sites in sample
        """          
        return self.thisptr.NumPoly() 
    def nummutations (self):
        """
        Number of mutations in sample.

        For this class, if there are k states at a site, there are k-1 mutations
        """          
        return self.thisptr.NumMutations() 
    def numsingletons (self):
        """
        Number of singletons (ancsetral or derived)
        """
        return self.thisptr.NumSingletons() 
    def numexternalmutations (self):
        """
        Number of derived mutations
        """
        return self.thisptr.NumExternalMutations() 
    #summary statistics of the site frequency spectrum
    def tajimasd (self):
        """
        Tajima's (1989) D statistics
        """
        return self.thisptr.TajimasD()
    def hprime (self,bint likeThorntonAndolfatto = False):
        """
        Normalized version of Fay and Wu's H statistic
        """
        return self.thisptr.Hprime(likeThorntonAndolfatto)
    def fulid (self):
        """
        Fu and Li's D
        """
        return self.thisptr.FuLiD()
    def fulif (self):
        """
        Fu and Li's F
        """
        return self.thisptr.FuLiF()
    def fulidstar (self):
        """
        Fu and Li's D*
        """
        return self.thisptr.FuLiDStar()
    def fulifstar (self):
        """
        Fu and Li's F*
        """
        return self.thisptr.FuLiFStar()
    def wallsb(self):
        """
        Jeff Wall's B statistic
        """
        return self.thisptr.WallsB()
    def wallsbprime(self):
        """
        Jeff Wall's B' statistic
        """
        return self.thisptr.WallsBprime()
    def wallsq(self):
        """
        Jeff Wall's Q statistic
        """
        return self.thisptr.WallsQ()
    def rm(self):
        """
        Hudson and Kaplan's lower bound on no. crossover events
        """
        return self.thisptr.Minrec()

cdef class polySIM:
    """
    "Factory" class to calculate summary statistics from binary data.

    This is a wrapper around libsequence's Sequence::PolySIM

    These objects are constructed from :class:`libsequence.polytable.simData`

    .. note:: 0/1 = ancestral/derived    
    """
    def __cinit__(self,simData d):
        self.thisptr = new PolySIM(dynamic_cast['SimData*'](d.thisptr))
    def __dealloc__(self):
        del self.thisptr
    def thetapi(self):
        """
        "Sum of site heterozygosity."  :math:`\\hat\\theta_\\pi = \\sum_i^S\\frac{c}{n}\\frac{n-c-1}{n-1}`,
        where :math:`S` is the number of polymorphisms and :math:`n` is the sample size.
        """
        return self.thisptr.ThetaPi()
    def thetaw (self):
        """
        Watterson's estimator of :math:`\\theta` from :math:`S`
        """
        return self.thisptr.ThetaW()
    def thetah (self):
        """
        Fay and Wu's estimator of :math:`\\theta`.
        """
        return self.thisptr.ThetaH()
    def thetal (self):
        """
        Normalized version of Fay and Wu's estimator of :math:`\\theta`.
        """
        return self.thisptr.ThetaL()

    #calculate various numbers related to polymorphism
    def numpoly(self):
        """
        Number of poylmorphic sites in sample
        """          
        return self.thisptr.NumPoly() 
    def nummutations (self):
        """
        Number of mutations in sample.

        For this class, this is the same as the number of polymorphic sites
        """          
        return self.thisptr.NumMutations() 
    def numsingletons (self):
        """
        Number of singletons (character states 0 or 1)
        """
        return self.thisptr.NumSingletons() 
    def numexternalmutations (self):
        """
        Number of derived mutations (character state 1)
        """
        return self.thisptr.NumExternalMutations() 
    #summary statistics of the site frequency spectrum
    def tajimasd (self):
        """
        Tajima's (1989) D statistics
        """
        return self.thisptr.TajimasD()
    def hprime (self,bint likeThorntonAndolfatto = False):
        """
        Normalized version of Fay and Wu's H statistic
        """
        return self.thisptr.Hprime(likeThorntonAndolfatto)
    def fulid (self):
        """
        Fu and Li's D
        """
        return self.thisptr.FuLiD()
    def fulif (self):
        """
        Fu and Li's F
        """
        return self.thisptr.FuLiF()
    def fulidstar (self):
        """
        Fu and Li's D*
        """
        return self.thisptr.FuLiDStar()
    def fulifstar (self):
        """
        Fu and Li's F*
        """
        return self.thisptr.FuLiFStar()
    def wallsb(self):
        """
        Jeff Wall's B statistic
        """
        return self.thisptr.WallsB()
    def wallsbprime(self):
        """
        Jeff Wall's B' statistic
        """
        return self.thisptr.WallsBprime()
    def wallsq(self):
        """
        Jeff Wall's Q statistic
        """
        return self.thisptr.WallsQ()
    def rm(self):
        """
        Hudson and Kaplan's lower bound on no. crossover events
        """
        return self.thisptr.Minrec()

##functions
def lhaf( polyTable pt, double l ):
    if isinstance(pt,simData):
        return lHaf(deref(dynamic_cast['SimData*'](pt.thisptr)),l)
    else:
        raise RuntimeError("lhaf: only simData objects are allowed")

def std_nSL(polyTable pt, double minfreq = 0., double binsize = 0.05, double[:] gmap = None):
    if isinstance(pt,simData):
        if gmap is None:
            return snSL(deref(dynamic_cast['SimData*'](pt.thisptr)),minfreq,binsize,NULL)
        else:
            return snSL(deref(dynamic_cast['SimData*'](pt.thisptr)),minfreq,binsize,&gmap[0])
    else:
        raise RuntimeError("lhaf: only simData objects are allowed")

def ld(polyTable p, bint haveOutgroup = False, unsigned outgroup = 0, unsigned mincount = 1,maxDist = None):
    """
    Pairwise "linkage disequilibrium" (LD) statistics

    :param p: A :class:`libsequence.polytable.polyTable`
    :param haveOutgroup: if True, then ougtroup is the index of the outroup sequence in p
    :param outgroup: The index of the outgroup sequence in p.  Not used if haveOutgroup is False
    :param mincount: Ignore site pairs where one or both sites have minor alleles occur < mincount times
    :param maxDist: Do not consider pairs of sites > maxDist units apart.  Default will be the maximum value of a C/C++ double.

    .. note:: This function skips sites with missing data, gaps, etc.
    """
    rv = {'i':[],
          'j':[],
          'rsq':[],
          'D':[],
          'Dprime':[]}
    if p.empty():
        return rv
    cdef unsigned i = 0
    cdef unsigned j = i + 1
    cdef vector[double] ldvals
    ldvals.resize(6)
    cdef numeric_limits[double] nl
    cdef double md = nl.max()
    if maxDist is not None:
        md = maxDist
    while Disequilibrium(p.thisptr,ldvals,&i,&j,haveOutgroup,outgroup,mincount,md):
        if ldvals[5] == 0:  ##site pair NOT skipped
            rv['i'].append(ldvals[0])
            rv['j'].append(ldvals[1])
            rv['rsq'].append(ldvals[2])
            rv['D'].append(ldvals[3])
            rv['Dprime'].append(ldvals[4])
    return rv
