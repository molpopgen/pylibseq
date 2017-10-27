Object lifetimes
===================================

This package provides two types of classes:

* Classes containing data, such as :class:`libsequence.polytable.SimData` and :class:`libsequence.polytable.PolySites`
* Classes analyzing the data in the previous type of class, *e.g.* :class:`libsequence.summstats.PolySIM` and :class:`libsequence.summstats.PolySNP`

In general, one uses these objects as follows:

.. code:: python

   #This is a 'data' class:
   data = SimData()
   data.assign( [(0.1,"010101"),(0.2,"101010")] )
   #This is an "analyze data" class:
   analyze = PolySim(data)
   d = analyze.tajimasd()

Objects that are 'analyze data' types **store pointers to the 'data' types**.  Thus, if you let the data go out of scope, delete it, reassign to it, or do any other sort of manipulation changing its internal state, then you will get a crash (segmentation fault) if you try to access members of the 'analyze data' types.
