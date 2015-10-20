pyseq: Python bindings for libsequence
***************************************************************

This package provides Python_ bindings for the C++11 library libsequence_.

The bindings are implemented using Cython_.

This package serves two roles:

* It provides a means of using some of the more widely-used bits of libsequence_ within the Python language
* The unit tests of this package also serve as unit tests for libsequence_.

What this package does **not** (currently) do:

* provide an interface for I/O operations.  Python I/O and C++ I/O are fundamentally very different.  Bridging the gap requires either adding features to Cython and/or adding modules to this package that depend on Boost.Python, which would add an additional C++ dependency to this package.

Requirements:
===================================

* libsequence_ must be installed on your system.  **Currently, this package requires the dev branch of libsequence**
* Python 2 or Python 3
* Cython_ must be installed on your system.  (Eventually, this will go away as an installation requirement).
* An up-to-date C++ compiler that is C++11 compatible via the flag -std=c++11.  Roughty, this means GCC >= 4.8 and clang >= 3.5.

You should install libsequence_ from source and Cython_ via your favorite Python package manager.

The supported platforms are Linux and OS X.

Installation:
=======================

.. code-block:: bash

   $ ./configure
   $ sudo python setup.py install

If you have libsequence in a "funny location" (*e.g.*, something other that /usr/local):

.. code-block:: bash

   $ CPPFLAGS=-I/path/to/libsequence/headers LDFLAGS=-L/path/to/libsequence/library sudo python setup.py install

For example, if libsequence is installed into /opt:

.. code-block:: bash

   $ CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib sudo python setup.py install

Unit testing:
=======================

.. code-block:: bash

   $ ./configure 
   $ python setup.py build_ext -i
   $ python -m unittest discover unit_test

Documentation:
======================

Eventually, documentation will be online.  Check the unit_test directory for working examples (well, working if the tests are passing!).

.. _libsequence: http://molpopgen.github.io/libsequence/
.. _Cython: http://www.cython.org/
.. _Python: http://www.cython.org/
