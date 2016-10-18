cdef class scheduler_init:
    def __cinit__(self,int nthreads = -1):
        self.init=unique_ptr[task_scheduler_init](new task_scheduler_init(nthreads))
    def terminate(self):
        if self.init.get().is_active():
            self.init.get().terminate()
    def change_nthreads(self,nthreads = None):
        if self.init.get().is_active():
            self.init.get().terminate()
        cdef int nt = self.init.get().default_num_threads()
        if nthreads is not None:
            nt=int(nthreads)
        self.init.get().initialize(nt)
