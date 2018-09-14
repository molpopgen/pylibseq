import msprime
import libsequence.variant_matrix as vmatrix
import libsequence.summstats as sstats
import seaborn as sns

D = []

for ts in msprime.simulate(100, num_replicates=100,
                           mutation_rate=25.0,
                           random_seed=42):
    m = vmatrix.VariantMatrix.from_TreeSequence(ts)
    ac = m.count_alleles()
    D.append(sstats.tajd(ac))

sns.distplot(D)
