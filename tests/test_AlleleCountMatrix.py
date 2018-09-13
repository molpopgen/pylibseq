import unittest
import msprime
import numpy as np
import libsequence.variant_matrix


class testAlleleCountMatrix(unittest.TestCase):
    @classmethod
    def setUp(self):
        ts = msprime.simulate(10, mutation_rate=10, random_seed=42)
        self.vm = libsequence.variant_matrix.VariantMatrix.from_TreeSequence(
            ts)
        self.ac = libsequence.variant_matrix.AlleleCountMatrix(self.vm)

    def testSubset(self):
        w = self.vm.window(0.25, 0.3)
        acw = libsequence.variant_matrix.AlleleCountMatrix(w)

        pa = np.array(self.vm.positions)
        p = np.where((pa >= 0.25) & (pa <= 0.3))[0]
        self.assertEqual(len(p), acw.nrow)

        ac_slice = self.ac[p.min():p.max()+1:1]
        self.assertTrue(np.array_equal(np.array(ac_slice),
                                       np.array(acw)))


if __name__ == "__main__":
    unittest.main()
