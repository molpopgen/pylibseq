Changes
===================================

Version 0.2.2
----------------------------------

* Package refactored into a single namespace
* Build system now uses cmake
* libsequence is included as a git submodule, removing all compilation dependencies except pybind11
* A VariantMatrix may now be constructed from numpy arrays without making copies
* Taking a window of a VariantMatrix is now copy-free

Version 0.2.1
----------------------------------

* Added :func:`libsequence.summstats.omega_max`
* Added :func:`libsequence.polytable.SimData.from_stdin`
