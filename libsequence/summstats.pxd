from libsequence.polytable cimport PolyTable,SimData,polyTable,simData
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport unique_ptr
from libcpp.unordered_map cimport unordered_map

cdef extern from "Sequence/PolySNP.hpp" namespace "Sequence" nogil:
    cdef cppclass PolySNP:
        PolySNP(const PolyTable * pt, const bool & haveOutgroup, const unsigned & outgroup, const bool & totMuts)

        #Summary stat fxns
        double ThetaPi() const                                #Tajima (1983)
        double ThetaW() const                            #Watterson's (1975) Theta
        double ThetaH() const                            #Theta from homozygosity (Fay and Wu 2001)
        double ThetaL() const                            #A variant on Fay and Wu's H
        #variances of estimators of 4Nu
        double VarPi() const
        double StochasticVarPi() const
        double SamplingVarPi() const
        double VarThetaW() const
        #calculate various numbers related to polymorphism
        unsigned NumPoly() const                                 #number of polymorphic sites in data
        unsigned NumMutations() const                    #number of inferred mutations 
        unsigned NumSingletons() const                   #number of mutants at frequency 1
        unsigned NumExternalMutations() const            #number of derived mutations 
        #summary statistics of the site frequency spectrum
        double TajimasD() const                          #Tajima's (1989) D
        double Hprime (const bool & likeThorntonAndolfatto) const  #A normalized statistic related to Fay and Wu's H
                      
        double Dnominator() const                        #Denominator of Tajima's D
        double FuLiD() const                             #Fu & Li's (1996) D
        double FuLiF() const                             #Fu & Li's (1996) F
        double FuLiDStar() const                         #Fu & Li's (1996) D*
        double FuLiFStar() const                         #Fu & Li's (1996) F*
        #summary statistics of haplotypes
        double DandVH() const                                    #Depaulis & Veuille (1998) Haplotype diversity
        unsigned DandVK() const                                  #Depaulis & Veuille (1998) number of haplotypes
        double WallsB() const                             #Jeff Wall's (1999) B  statistic
        unsigned WallsBprime() const                      #Jeff Wall's (1999) B' statistic
        double WallsQ() const                             #Jeff Wall's (1999) Q  statistic
        #recombination
        double HudsonsC() const                                  #Dick Hudson's (1987) Chat = 4Nr
        unsigned Minrec() const    

cdef extern from "Sequence/PolySIM.hpp" namespace "Sequence" nogil:
    cdef cppclass PolySIM:
        PolySIM(const SimData * pt)

        #Summary stat fxns
        double ThetaPi() const
        double ThetaW () const
        double ThetaH () const
        double ThetaL () const

        #calculate various numbers related to polymorphism
        unsigned NumPoly() const                                 
        unsigned NumMutations () const
        unsigned NumSingletons () const
        unsigned NumExternalMutations () const
        #summary statistics of the site frequency spectrum
        double TajimasD () const
        double Hprime (const bool & likeThorntonAndolfatto) const
        double Dnominator () const
        double FuLiD () const
        double FuLiF () const
        double FuLiDStar () const
        double FuLiFStar () const
        double WallsB() const
        unsigned WallsBprime() const
        double WallsQ() const
        #Hudson's Haplotype Partition Test
        int HudsonsHaplotypeTest (const int & subsize,const int & subss) const
        #recombination
        unsigned Minrec () const

cdef class polySNP:
    cdef unique_ptr[PolySNP] thisptr

cdef class polySIM:
    cdef unique_ptr[PolySIM] thisptr

##Functions from libseq
cdef extern from "Sequence/SummStats/Garud.hpp" namespace "Sequence" nogil:
    cdef cppclass GarudStats:
        double H1
        double H12
        double H2H1

    GarudStats H1H12(const SimData &)

cdef extern from "Sequence/SummStats/lHaf.hpp" namespace "Sequence" nogil:
    vector[double] lHaf( const SimData & data, const double l )

cdef extern from "Sequence/SummStats/nSL.hpp" namespace "Sequence" nogil:
    #These functions can throw exceptions when maps are used that do not
    #contain positions in d.
    vector[pair[double,double]] nSL_t(const SimData & d, const unordered_map[double,double] & gmap) except +
    pair[double,double] snSL(const SimData & d,const double minfreq, const double binsize, const unordered_map[double,double] & gmap) except +

cdef extern from "Sequence/Recombination.hpp" namespace "Sequence" nogil:
    cdef struct PairwiseLDstats:
        double i
        double j
        double rsq
        double D
        double Dprime
        bint skipped

cdef extern from "Sequence/Recombination.hpp" namespace "Sequence::Recombination" nogil:
    vector[PairwiseLDstats] Disequilibrium(const PolyTable *,const bint & haveOutgroup,const unsigned & outgroup,const unsigned & mincount,
            const double max_distance)
