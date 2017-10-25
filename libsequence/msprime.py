from libsequence.polytable import SimData


def make_SimData(g):
    """
    Construct a :class:`libsequence.polytable.SimData` from
    the output of msprime.

    :param g: The output from msprime

    .. note::
        Thanks to Jerome Kelleher for pointing out the
        quick implementation using msprime >= 0.4.0.

    .. versionchanged:: 0.2.0
        Much faster implementation via haplotypes in tree sequence

    Example:

    >>> import msprime as msp
    >>> from libsequence.msprime import make_SimData
    >>> g = msp.simulate(sample_size = 10, Ne=1e6,\
                         recombination_rate=1e-8, \
                         mutation_rate=1e-8,length=1e4)
    >>> s = make_SimData(g)
    """
    return SimData([i.position for i in g.mutations()],
                   [i for i in g.haplotypes()])
