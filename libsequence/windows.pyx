from libsequence.polytable cimport polyTable,polySites,simData
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libcpp.cast cimport dynamic_cast
from cython.operator cimport dereference as deref

cdef fill_from_SimData(const SimData *d ,double window_size, double step_len, double starting_pos, double ending_pos):
    cdef PolyTableSlice[SimData] * windows=new PolyTableSlice[SimData](d.sbegin(),d.send(),window_size,step_len,starting_pos,ending_pos)
    cdef SimData d2
    cdef vector[pair[double,string]] temp
    wins = []
    for i in range(windows.size()):
        d2=deref(windows)[i]
        temp.assign(d2.sbegin(),d2.send())
        wins.append(simData(temp))
    del windows
    return wins

cdef fill_from_PolySites(const PolySites *d ,double window_size, double step_len, double starting_pos, double ending_pos):
    cdef PolyTableSlice[PolySites] * windows=new PolyTableSlice[PolySites](d.sbegin(),d.send(),window_size,step_len,starting_pos,ending_pos)
    cdef PolySites d2
    cdef vector[pair[double,string]] temp
    wins = []
    for i in range(windows.size()):
        d2=deref(windows)[i]
        temp.assign(d2.sbegin(),d2.send())
        wins.append(simData(temp))
    del windows
    return wins
    
cdef class Windows:
    """
    An iterable list of sliding windows created from a :class:`libsequence.polytable.polyTable`
    """
    def __cinit__(self, polyTable pt, double window_size, double step_len, double starting_pos = 0., double ending_pos = 1):
        if isinstance(pt,simData):
            self.windows = fill_from_SimData(dynamic_cast['SimData*'](pt.thisptr.get()),window_size,step_len,starting_pos,ending_pos)
        else:
            self.windows = fill_from_PolySites(dynamic_cast['PolySites*'](pt.thisptr.get()),window_size,step_len,starting_pos,ending_pos)
    def __dealloc__(self):
        self.windows=[]
    def __iter__(self):
        return iter(self.windows)
    def __next__(self):
        return next(self.windows)
    def __getitem__(self,i):
        return self.windows[i]
    def __len__(self):
        return len(self.windows)
    
    
