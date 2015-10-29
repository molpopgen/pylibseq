from libsequence.polytable cimport polySites,simData
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from cython.operator cimport dereference as deref

cdef class simDataWindows:
    """
    Calculate sliding windows from a :class:`libsequence.polytable.simData`
    """
    def __cinit__(self, simData d, double window_size, double step_len, double starting_pos = 0., double ending_pos = 1):
        self.windows = new PolyTableSlice[SimData](d.thisptr.sbegin(),
                                                   d.thisptr.send(),
                                                   window_size,step_len,starting_pos,ending_pos)
        self.wins = list()
        cdef SimData d2
        cdef vector[pair[double,string]] temp
        for i in range(self.windows.size()):
            d2 = deref(self.windows)[i]
            temp.assign(d2.sbegin(),d2.send())
            self.wins.append(simData(temp))
    def __dealloc__(self):
        del self.windows
        self.wins=[]
    def __iter__(self):
        return iter(self.wins)
    def __next__(self):
        return next(self.wins)
    def __getitem__(self,i):
        return self.wins[i]
    def __len__(self):
        return self.windows.size()

cdef class polySitesWindows:
    """
    Calculate sliding windows from a :class:`libsequence.polytable.polySites`
    """
    def __cinit__(self, polySites d, double window_size, double step_len, double starting_pos = 0., double ending_pos = 1):
        self.windows = new PolyTableSlice[PolySites](d.thisptr.sbegin(),
                                                     d.thisptr.send(),
                                                     window_size,step_len,starting_pos,ending_pos)
        self.wins = list()
        cdef PolySites d2
        cdef vector[pair[double,string]] temp
        for i in range(self.windows.size()):
            d2 = deref(self.windows)[i]
            temp.assign(d2.sbegin(),d2.send())
            self.wins.append(polySites(temp))
    def __dealloc__(self):
        del self.windows
        self.wins=[]
    def __iter__(self):
        return iter(self.wins)
    def __next__(self):
        return next(self.wins)
    def __getitem__(self,i):
        return self.wins[i]
    def __len__(self):
        return self.windows.size()
    
