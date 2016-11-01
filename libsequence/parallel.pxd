from libcpp.memory cimport unique_ptr

cdef extern from "tbb/task_scheduler_init.h" namespace "tbb" nogil:
    cdef cppclass task_scheduler_init:
        #this constructor wrapper is ignoring
        #the second argument, which defaults to 0
        task_scheduler_init(int max_threads)
        void initialize(int max_threads)
        void terminate()
        @staticmethod
        int default_num_threads()
        bint is_active() const

cdef extern from "Sequence/threading.hpp" namespace "Sequence" nogil:
    unique_ptr[task_scheduler_init] init_tbb(int nthreads) except +

cdef class scheduler_init(object):
    cdef unique_ptr[task_scheduler_init] init

