.. _summstatsVM:

Calculating summary statistics from a VariantMatrix
================================================================

Summary statistic functions are found in :mod:`libsequence.summstats`.  Certain 
summary statistics can be calculated simply from allele counts via
:class:`libsequence.variant_matrix.AlleleCountMatrix`. Calculations invovling linkage disequilibrium
will typically require the entire :class:`libsequence.variant_matrix.VariantMatrix` object.

Distribution of Tajima's D from msprime
------------------------------------------------------------------------------

.. plot:: pyplots/tajd.py
    :include-source:


Binning the :math:`nS_L` statistic from coalescent simulations
------------------------------------------------------------------------------

"Haplotype homozygosity" statistics are very popular in certain circles.  Doing any kind of statistics
with them requires binning "neutral" or "reference" SNPs by allele frequency so that you can convert
raw scores into z-scores.  In simulation studies, it makes sense to calibrate the means and standard
deviations per bin based on simulations of a null model.  Here, I show how to do this for the :math:`nS_L`
statistic of :cite:`Ferrer-Admetlla2014-wa`.

The relevant function here is :func:`libsequence.summstats.nsl`, which returns instances of
:class:`libsequence.summstats.nSLresults`.


.. plot:: pyplots/nSLbins.py
    :include-source:
