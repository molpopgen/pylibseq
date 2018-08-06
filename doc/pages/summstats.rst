.. _summarrystats:

Calculating summary statistics
================================================================

This page is a tutorial on calculating summary statistics from a :class:`libsequence.variant_matrix.VariantMatrix`.  We
also show how to interchange data between msprime_ :cite:`Kelleher2016-cb` and pylibseq.

.. ipython:: python

    import msprime
    import libsequence.variant_matrix as vmat
    import libsequence.summstats as sstats
    import numpy as np

    # Get a TreeSequence from msprime
    ts = msprime.simulate(100, mutation_rate=250., recombination_rate=250., random_seed=42)

    # TreeSequence -> VariantMatrix
    vm = vmat.VariantMatrix.from_TreeSequence(ts)

    # Standard summary stats
    pi = sstats.thetapi(vm)
    tajd = sstats.tajd(vm)
    hprime = sstats.hprime(vm, 0) # 0 = reference state = ancestral state

Some more complex descriptors of the data are available

.. ipython:: python

    diffs = sstats.difference_matrix(vm)

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

    diff_yes_or_no = sstats.is_different_matrix(vm)

The contents of this matrix have the exact same layout as `diffs` described above.  The difference is that the data
elements are encoded as 0 = identical, 1 = different.  This calculation is **much** faster than the previous.

.. note::

    Missing data do not contribute to samples being considered different, nor to the 
    number of differences.

It is also possible to get a unique label assigned to each haplotype:

.. ipython:: python

    labels = np.array(sstats.label_haplotypes(vm),dtype=np.int32)
    print(len(np.unique(labels)))

Internally, the results from `is_different_matrix` are used to assign the labels.

These labels are likewise used internally to count the number of haplotypes:

.. ipython:: python

    print(sstats.number_of_haplotypes(vm))
    # Confirm result via direct comparison to 
    # the data from msprime:
    print(len(np.unique(gm.transpose(),axis=0)))

What about performance?

.. ipython:: python

    %timeit -n 10 -r 10 sstats.number_of_haplotypes(vm)

.. ipython:: python
    
    %timeit -n 10 -r 10 len(np.unique(gm.transpose(),axis=0))
   

.. _msprime: http://msprime.readthedocs.io

