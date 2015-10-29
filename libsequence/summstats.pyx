from libcpp.cast cimport dynamic_cast
from cython.view cimport array as cvarray
from cpython cimport array
from cython.operator cimport dereference as deref

cdef extern from "<limits>" namespace "std" nogil:
    cdef cppclass numeric_limits[T]:
        @staticmethod
        T max()
        
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
    """
    :math:`l-HAF` from Ronen et al. DOI:10.1371/journal.pgen.1005527
    
    :param pt: A :class:`libsequence.polytable.polyTable`
    
    :return: The :math:`l-HAF` statistic for each haplotype in pt
    
    :rtype: list
    
    .. note:: Only :class:`libsequence.polytable.simData` types currently supported
    """
    if isinstance(pt,simData):
        return lHaf(deref(dynamic_cast['SimData*'](pt.thisptr)),l)
    else:
        raise RuntimeError("lhaf: only simData objects are allowed")

def nSLiHS(polyTable pt, double[:] gmap = None):
    """
    "Raw"/unstandardizes :math:`nS_L` and iHS from Ferrer-Admetlla et al. doi:10.1093/molbev/msu077.

    :param pt: A :class:`libsequence.polytable.polyTable`

    :return: A list of (nSL,iHS) tuples

    :rtype: list
    
    .. note:: Only :class:`libsequence.polytable.simData` types currently supported
    """
    if isinstance(pt,simData):
        rv = []
        if gmap is None:
            for core in range(pt.numsites()):
                rv.append(nSL(core,deref(dynamic_cast['SimData*'](pt.thisptr)),NULL))
        else:
            for core in range(pt.numsites()):
                rv.append(nSL(core,deref(dynamic_cast['SimData*'](pt.thisptr)),&gmap[0]))
        return rv
    else:
        raise RuntimeError("nSL: only simData objects are allowed")
    
def std_nSLiHS(polyTable pt, double minfreq = 0., double binsize = 0.05, double[:] gmap = None):
    """
    Standardized :math:`nS_L` statistic from Ferrer-Admetlla et al. doi:10.1093/molbev/msu077

    :param pt: A :class:`libsequence.polytable.polyTable`
    :param minfreq: Ignore markers with frequency < this value
    :param binsize: Standardize statistic in frequency bings of this width
    :parma gmap: A genetic map.  The length of this array should be equal to the number of sites in pt, and represent the positions of those sites on the genetic map.

    :return: A tuple. The first value is max standardized nSL over all bins.  The second is max iHS over all bins, where iHS is calculated according to Ferrer-Admetlla et al.
    The maxmimums are calculated based on absolute value.

    :rtype: float

    .. note:: Only :class:`libsequence.polytable.simData` types currently supported
    """
    if isinstance(pt,simData):
        if gmap is None:
            return snSL(deref(dynamic_cast['SimData*'](pt.thisptr)),minfreq,binsize,NULL)
        else:
            return snSL(deref(dynamic_cast['SimData*'](pt.thisptr)),minfreq,binsize,&gmap[0])
    else:
        raise RuntimeError("std_nSL: only simData objects are allowed")

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
    cdef double md = numeric_limits[double].max()
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

def garudStats(polyTable pt):
   """
   Statistics from Garud et al. doi:10.1371/journal.pgen.1005004.g011

   :param pt: A :class:`libsequence.polytable.polyTable`

   :return: The H1, H12, and H2/H1 statistics

   :rtype: dict

   .. note:: Only :class:`libsequence.polytable.simData` types currently supported
   """
   cdef GarudStats stats
   if isinstance(pt,simData):
       stats = H1H12(deref(dynamic_cast['SimData*'](pt.thisptr)))
       return {"H1":stats.H1,"H12":stats.H12,"H2H1":stats.H2H1}
   else:
       raise RuntimeError("garudStats: only simData objects are allowed")       

