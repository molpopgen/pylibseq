import unittest
import libsequence
import cpptests
import numpy as np


class TestAlleleCountMatrix(unittest.TestCase):
    def test_create(self):
        am = cpptests.create_AlleleCountMatrix()
        self.assertTrue(isinstance(am, libsequence.AlleleCountMatrix))


class TestVariantMatrix(unittest.TestCase):
    def test_create(self):
        data = [1, 0, 1, 0]
        pos = [0.1, 0.2]
        vm = cpptests.create_VariantMatrix(data, pos)
        self.assertTrue(np.array_equal(np.array(data).reshape(2, 2),
                                       np.array(vm.data)))
        self.assertTrue(np.array_equal(np.array(pos), np.array(vm.positions)))


if __name__ == "__main__":
    unittest.main()
