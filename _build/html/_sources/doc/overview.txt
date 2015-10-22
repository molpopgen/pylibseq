
Example workflow
================

We'll work with pyseq's wrapper to libsequence's SimData, which is used
to process bi-allele data encoded as 0/1 = ancestral/derived,
respectively

.. code:: python

    from __future__ import print_function
    from libsequence.polytable import simData

Assigning to an object
----------------------

You may assign from a list of tuples.

Each tuple is a site: (pos:genotypes)

Here, there are 2 sites and a sample size of :math:`n=4`

.. code:: python

    rawData1 = [(0.1,'0101'),(0.2,'1010')]

.. code:: python

    sd = simData()

.. code:: python

    sd.assign(rawData1)
    sd.numsites()




.. parsed-literal::

    2



.. code:: python

    sd.size()




.. parsed-literal::

    4



.. code:: python

    sd.pos()




.. parsed-literal::

    [0.1, 0.2]



.. code:: python

    sd.data()




.. parsed-literal::

    ['01', '10', '01', '10']



Or, you can assign from separate list of positions and haplotypes

.. code:: python

    rawDataPos = [0.1,0.2]
    rawDataGenos = ['01','10','01','10']
    sd.assign_sep(rawDataPos,rawDataGenos)

.. code:: python

    sd.numsites()




.. parsed-literal::

    2



.. code:: python

    sd.size()




.. parsed-literal::

    4



.. code:: python

    sd.pos()




.. parsed-literal::

    [0.1, 0.2]



.. code:: python

    sd.data()




.. parsed-literal::

    ['01', '10', '01', '10']



Summary statistics
------------------

Let's calculate some basic summary statistics

See :class:`libsequence.summstats.polySIM` for more documentation

.. code:: python

    from libsequence.summstats import polySIM
    #ms 10 1 -s 10 -I 2 5 5 0.05
    rawDataPos=[0.0997, 0.2551, 0.3600, 0.4831, 0.5205, 0.5668, 0.5824, 0.6213, 0.7499, 0.9669]
    rawDataGenos=['0000001010',
                  '0000000011',
                  '0000001010',
                  '0000001010',
                  '0000001010',
                  '1111010100',
                  '1111010100',
                  '1111110100',
                  '1111010100',
                  '1111010100']
    sd.assign_sep(rawDataPos,rawDataGenos)

.. code:: python

    ps = polySIM(sd)

.. code:: python

    ps.thetapi()




.. parsed-literal::

    4.822222222222222



.. code:: python

    ps.thetaw()




.. parsed-literal::

    3.5348576237901534



.. code:: python

    ps.tajimasd()




.. parsed-literal::

    1.6142469967484658



Sliding windows
---------------

.. code:: python

    from libsequence.windows import simDataWindows

.. code:: python

    w = simDataWindows(sd,window_size=0.1,step_len=0.05,starting_pos=0.,ending_pos=1.0)

.. code:: python

    len(w)




.. parsed-literal::

    20



.. code:: python

    for i in range(len(w)):
        #Each window is a simData
        wi = w[i]
        pswi = polySIM(wi)
        print(pswi.thetaw())


.. parsed-literal::

    0.353485762379
    0.353485762379
    0.0
    0.0
    0.353485762379
    0.353485762379
    0.353485762379
    0.353485762379
    0.353485762379
    0.706971524758
    1.06045728714
    1.06045728714
    0.353485762379
    0.353485762379
    0.353485762379
    0.0
    0.0
    0.0
    0.353485762379
    0.353485762379


:math:`F_{ST}`
--------------

Let's pretend that our data are from two demes of sizes n/2 each.

Note that most flavors of :math:`F_{ST}` are very similar to one
another. See Charlesworth, B. (1998) Mol. Biol. Evol. 15(5): 538-543 for
a great overview.

.. code:: python

    from libsequence.fst import fst
    sd.size()
    f = fst(sd,[5,5])

.. code:: python

    #Hudson, Slatkin, and Maddison's FST:
    f.hsm()




.. parsed-literal::

    0.9268292682926829



.. code:: python

    #Slatkin's
    f.slatkin()




.. parsed-literal::

    0.8636363636363636



.. code:: python

    #Hudson, Boos, and Kaplan, which is also Nei's Gst:
    f.hbk()




.. parsed-literal::

    0.8636363636363635



.. code:: python

    #Positions of snps shared b/w demes 0 and 1
    f.shared(0,1)




.. parsed-literal::

    []



.. code:: python

    #Positions of private mutations in deme 0 and 1:
    f.priv(0,1)




.. parsed-literal::

    {0: [0.5824, 0.9669], 1: [0.5205]}



.. code:: python

    #Positions of fixed differences between demes 0 and 1:
    f.fixed(0,1)




.. parsed-literal::

    [0.0997, 0.2551, 0.36, 0.4831, 0.5668, 0.6213, 0.7499]


