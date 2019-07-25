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

    def test_create_from_external_cpp_vectors(self):
        vm = cpptests.create_VariantMatrix_from_cpp_vectors()
        edata = np.array([0, 1, 1, 0], dtype=np.int8).reshape(2, 2)
        epos = np.array([0.1, 0.2])
        self.assertTrue(np.array_equal(vm.data, edata))
        self.assertTrue(np.array_equal(vm.positions, epos))


if __name__ == "__main__":
    unittest.main()
