Parallel computation
================================

Under the hood, libsequence using Intel's TBB_ library for parallelism.  By default, TBB_ will automatically deterimine the optimal number of cores to use.  However, that behavior may not be desirable in certain contexts.  For example, on managed HPC systems you may be expected to keep your use of cores under control.  You make limit the maximum number of threads used as follows:

>>> from libsequence.parallel import Scheduler
>>> init = Scheduler(10)

Any value less than one will result in the TBB_ default number of threads being used.

.. _TBB: http://www.threadbuildingblocks.org
