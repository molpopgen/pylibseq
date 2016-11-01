cdef class SimDataVec:
    """
    Wrapper for std::vector<Sequence::SimData> (not the Cython extension 
    class libsequence.polytable.simData!!)

    This type is useful in the event that you want to process data
    sets in parallel using OpenMP.    

    .. note:: The use case for this type is a Cython-based package extending pylibseq.
    """
    def __cinit__(self):
        self.vec = vector[SimData]()
    def __dealloc__(self):
        self.vec.clear()
    def __init__(self,const vector[polySiteVector] & p):
        """
        Constructor.

        :param p: A list of a list of tuples representing the variation data.

        Example:

        >>> import libsequence.extensions as le
        >>> #Create a vector of 5 identical data sets
        >>> #with 2 variable sites per data set.
        >>> x = [[(0.1,"00100"),(0.2,"11000")]]*5
        >>> y = le.SimDataVec(x)
        """
        cdef int i=0
        cdef int n = p.size()
        cdef psite_vec_itr a,b
        while i<n:
            a=p[i].begin()
            b=p[i].end()
            self.vec.push_back(SimData(a,b))
            i+=1
