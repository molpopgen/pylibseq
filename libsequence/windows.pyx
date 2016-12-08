from libsequence.polytable cimport PolyTable,CppPolySites,CppSimData
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.cast cimport dynamic_cast
from libcpp.memory cimport unique_ptr 
from cython.operator cimport dereference as deref

cdef fill_from_SimData(const CppSimData *d ,double window_size, double step_len, double starting_pos, double ending_pos):
    cdef unique_ptr[PolyTableSlice[CppSimData]] windows
    windows.reset(new PolyTableSlice[CppSimData](d.sbegin(),d.send(),window_size,step_len,starting_pos,ending_pos))
    cdef CppSimData d2
    cdef vector[pair[double,string]] temp
    wins = []
    for i in range(windows.get().size()):
        d2=deref(windows)[i]
        temp.assign(d2.sbegin(),d2.send())
        wins.append(SimData(temp))
    return wins

cdef fill_from_PolySites(const CppPolySites *d ,double window_size, double step_len, double starting_pos, double ending_pos):
    cdef unique_ptr[PolyTableSlice[CppPolySites]] windows
    windows.reset(new PolyTableSlice[CppPolySites](d.sbegin(),d.send(),window_size,step_len,starting_pos,ending_pos))
    cdef CppPolySites d2
    cdef vector[pair[double,string]] temp
    wins = []
    for i in range(windows.get().size()):
        d2=deref(windows)[i]
        temp.assign(d2.sbegin(),d2.send())
        wins.append(PolySites(temp))
    return wins
    
cdef class Windows:
    """
    An iterable list of sliding windows created from a :class:`libsequence.polytable.polyTable`
    """
    def __cinit__(self, PolyTable pt, double window_size, double step_len, double starting_pos = 0., double ending_pos = 1):
        if isinstance(pt,SimData):
            self.windows = fill_from_SimData(dynamic_cast['CppSimData*'](pt.thisptr.get()),window_size,step_len,starting_pos,ending_pos)
        else:
            self.windows = fill_from_PolySites(dynamic_cast['CppPolySites*'](pt.thisptr.get()),window_size,step_len,starting_pos,ending_pos)
    def __iter__(self):
        return iter(self.windows)
    def __next__(self):
        return next(self.windows)
    def __getitem__(self,i):
        return self.windows[i]
    def __len__(self):
        return len(self.windows)
    
    
