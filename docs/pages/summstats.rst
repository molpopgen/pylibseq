.. _summarrystats:

Calculating summary statistics
================================================================

This page is a tutorial on calculating summary statistics from a :class:`libsequence.VariantMatrix`.  We
also show how to interchange data between msprime_ :cite:`Kelleher2016-cb` and pylibseq.

Simple overview of concepts
-------------------------------------------------------

.. ipython:: python

    import msprime
    import libsequence
    import numpy as np

    # Get a TreeSequence from msprime
    ts = msprime.simulate(100, mutation_rate=250., recombination_rate=250., random_seed=42)

    # TreeSequence -> VariantMatrix
    vm = libsequence.VariantMatrix.from_TreeSequence(ts)
    # VariantMatrix -> AlleleCountMatrix
    ac = vm.count_alleles()

Certain summary statistics can be calculated simply from allele counts via
:class:`libsequence.AlleleCountMatrix`. Calculations invovling linkage disequilibrium
will typically require the entire :class:`libsequence.VariantMatrix` object.

The following examples rely on msprime :cite:`Kelleher2016-cb` to generate input data, but 
the concepts hold more generally.

List of "classic"/standard summary statistics
-------------------------------------------------------

Nucleotide diversity, or :math:`\pi` :cite:`Tajima1983-it`:

.. autofunction:: libsequence.thetapi

.. ipython:: python

    print(libsequence.thetapi(ac))

Watterson's estimator of :math:`\theta`, :cite:`Watterson1975-ej`:

.. autofunction:: libsequence.thetaw

.. ipython:: python

    print(libsequence.thetaw(ac))

Tajima's :math:`D`, :cite:`Tajima1989-de`:

.. autofunction:: libsequence.tajd

.. ipython:: python

    print(libsequence.tajd(ac))

Fay and Wu's test statistic, :math:`\pi - \hat\theta_H` :cite:`Fay2000-ef` has two overloads.
The first takes
a single ancestral state as an argument.  The latter takes a
list of ancestral states (one per site in the variant matrix):

.. autofunction:: libsequence.faywuh

.. ipython:: python

    print(libsequence.faywuh(ac, 0))

The :math:`H'` statistic :cite:`Zeng2006-is` should be preferred to
the Fay and Wu test statistic:

.. autofunction:: libsequence.hprime

.. ipython:: python

    print(libsequence.hprime(ac, 0))

Summaries of haplotype variation :cite:`Depaulis1998-ol`:

.. autofunction:: libsequence.number_of_haplotypes
.. autofunction:: libsequence.haplotype_diversity

.. ipython:: python

    print(libsequence.number_of_haplotypes(vm))
    print(libsequence.haplotype_diversity(vm))

The :math:`nS_L` and :math:`iHS` "haplotype homozygosity" statistics
defined in :cite:`Ferrer-Admetlla2014-wa`:

.. autofunction:: libsequence.nsl

This function returns a list of instances of the following:

.. autoclass:: libsequence.nSLresults
    :members:

Values of `nan` imply that haplotype homozygosity for a particular core
variant was not broken within the matrix.  This is the same procedure
used in :cite:`Ferrer-Admetlla2014-wa`.  The return value order
corresponds to the row order of the input data:

.. ipython:: python

    nsl = libsequence.nsl(vm, 0)
    for i in nsl[:4]:
       print(i.core_count, i.ihs, i.nsl)

The following variant only allows haplotype homozygosity to be
broken up by mutations where the non-reference count is :math:`\leq x`:

.. autofunction:: libsequence.nslx

.. ipython:: python

    nslx = libsequence.nslx(vm, 0, 1)

The haplotype diversity statistics from :cite:`Garud2015-ob`:

.. autofunction:: libsequence.garud_statistics

.. autoclass:: libsequence.GarudStats
    :members:

.. ipython:: python

    g = libsequence.garud_statistics(vm)
    print(g.H1, g.H12, g.H2H1)


.. autofunction:: libsequence.two_locus_haplotype_counts
.. autofunction:: libsequence.allele_counts
.. autofunction:: libsequence.non_reference_allele_counts
.. autoclass:: libsequence.AlleleCounts
    :members:

Distribution of Tajima's D from msprime 
------------------------------------------------------------------------------

.. plot:: pyplots/tajd.py
    :include-source:

"Sliding window" analysis
----------------------------------------------------------------------------

To analyze different genomic intervals, or "windows", you have two options.
Any statistic needing the allele count data can easily be gotten via a slice 
of your existing data:

.. ipython:: python

    lefts = np.arange(0, 1, 0.2)
    for l in lefts:
        p = np.where((vm.positions >= l)&(vm.positions < l+0.2))[0]
        ac_slice=ac[p[0]:p[-1]+1]
        assert len(ac_slice) == len(p)
        print(libsequence.thetapi(ac_slice))

We may also use :func:`libsequence.VariantMatrix.window` to obtain new `VariantMatrix` instances
corresponding to sub-windows of the data:

.. ipython:: python

    for l in lefts:
        w = vm.window(l, l+0.2)
        acw = w.count_alleles()
        print(libsequence.thetapi(acw))

Window creation is :math:`O(log(vm.nsites))` in time and has trivial additional memory requirements,
as the returned object does not own its own data buffer.

Other useful statistics
----------------------------------------------------------------

Some more complex descriptors of the data are available

.. ipython:: python

    diffs = libsequence.difference_matrix(vm)

Here, `diffs` contains data for representing the distance matrix for the data.  Specifically, these values can be used
to fill the upper triangle of a matrix:

.. ipython:: python

    dm = np.array(diffs, dtype=np.int32)
    m = np.array([0]*(vm.nsam*vm.nsam))
    # Get the indices of the upper-diagonal
    # of an nsam*nsam matrix, ignoring
    # the diagonal:
    idx = np.triu_indices(vm.nsam,k=1)
    m=m.reshape((vm.nsam,vm.nsam))
    m[idx]=dm
    print(m)

    # Let's confirm our results via brute-force:
    gm = ts.genotype_matrix()
    dummy = 0
    for i in range(gm.shape[1]-1):
        for j in range(i+1,gm.shape[1]):
            diffs = np.where(gm[:,i] != gm[:,j])
            assert(len(diffs[0])==dm[dummy])
            dummy+=1

In a similar fashion, we can obtain true/false data on whether pairs of haplotypes differ:

.. ipython:: python

    diff_yes_or_no = libsequence.is_different_matrix(vm)

.. autofunction:: libsequence.is_different_matrix

The contents of this matrix have the exact same layout as `diffs` described above.  The difference is that the data
elements are encoded as 0 = identical, 1 = different.  This calculation is **much** faster than the previous.

.. note::

    Missing data do not contribute to samples being considered different, nor to the 
    number of differences.

It is also possible to get a unique label assigned to each haplotype:

.. ipython:: python

    labels = np.array(libsequence.label_haplotypes(vm),dtype=np.int32)
    print(len(np.unique(labels)))

.. autofunction:: libsequence.label_haplotypes

Internally, the results from `is_different_matrix` are used to assign the labels.

These labels are likewise used internally to count the number of haplotypes:

.. ipython:: python

    print(libsequence.number_of_haplotypes(vm))
    # Confirm result via direct comparison to 
    # the data from msprime:
    print(len(np.unique(gm.transpose(),axis=0)))

What about performance?

.. ipython:: python

    %timeit -n 10 -r 10 libsequence.number_of_haplotypes(vm)

.. ipython:: python
    
    %timeit -n 10 -r 10 len(np.unique(gm.transpose(),axis=0))
   

.. _msprime: http://msprime.readthedocs.io

