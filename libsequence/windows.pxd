from libsequence.polysitevector cimport psite_vec_const_itr
from libsequence.polytable cimport PolySites,SimData

cdef extern from "Sequence/PolyTableSlice.hpp" namespace "Sequence":
    cdef cppclass PolyTableSlice[T]:
        PolyTableSlice(psite_vec_const_itr beg,psite_vec_const_itr end,
                       const double & window_size,
                       const double & step_len,
                       const double & starting_pos,
                       const double & ending_pos) except +

        T operator[](const unsigned &) const
        unsigned size() const

ctypedef PolyTableSlice[SimData] SimDataWindows
ctypedef PolyTableSlice[PolySites] PolySitesWindows

cdef class Windows:
    cdef public object windows
        
