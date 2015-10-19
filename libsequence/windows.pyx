from libsequence.polytable cimport polySites,simData
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from cython.operator cimport dereference as deref

cdef class simDataWindows:
    def __cinit__(self, simData d, double window_size, double step_len, double starting_pos = 0., double ending_pos = 1):
        self.windows = new PolyTableSlice[SimData](d.thisptr.sbegin(),
                                                   d.thisptr.send(),
                                                   window_size,step_len,starting_pos,ending_pos)
    def __dealloc__(self):
        del self.windows
    def __getitem__(self,i):
        cdef SimData d = deref(self.windows)[i]
        cdef vector[pair[double,string]] temp;
        temp.assign(d.sbegin(),d.send())
        rv = simData()
        rv.assign(temp)
        return rv
    def __len__(self):
        return self.windows.size()

cdef class polySitesWindows:
    def __cinit__(self, polySites d, double window_size, double step_len, double starting_pos = 0., double ending_pos = 1):
        self.windows = new PolyTableSlice[PolySites](d.thisptr.sbegin(),
                                                     d.thisptr.send(),
                                                     window_size,step_len,starting_pos,ending_pos)
    def __dealloc__(self):
        del self.windows
    def __getitem__(self,i):
        cdef PolySites d = deref(self.windows)[i]
        cdef vector[pair[double,string]] temp;
        temp.assign(d.sbegin(),d.send())
        rv = simData()
        rv.assign(temp)
        return rv
    def __len__(self):
        return self.windows.size()
    
