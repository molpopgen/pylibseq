from libsequence.polytable import simData
def make_simData(g):
    """
    Construct a :class:`libsequence.polytable.simData` from 
    the output of msprime.

    :param g: The output from msprime

    .. note:: Thanks to Jerome Kelleher for pointing out the quick implementation using msprime >= 0.4.0.

    Example:

    >>> import msprime as msp
    >>> from libsequence.msprime import make_simData
    >>> g = msp.simulate(sample_size = 10,Ne=1e6, recombination_rate=1e-8,mutation_rate=1e-8,length=1e4)
    >>> s = make_simData(g)
    """
    return simData([(v.position, v.genotypes) for v in g.variants(as_bytes=True)])
