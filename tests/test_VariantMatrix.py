import unittest
import libsequence.variant_matrix
import numpy as np

V8 = libsequence.variant_matrix.VectorInt8
VD = libsequence.variant_matrix.VectorDouble


class testVariantMatrix(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.dlist = [0, 1, 1, 0, 0, 0, 0, 1]
        self.plist = [0.1, 0.2]
        self.data = V8(self.dlist)
        self.pos = VD(self.plist)
        self.m = libsequence.variant_matrix.VariantMatrix(self.data, self.pos)
        self.m = libsequence.variant_matrix.VariantMatrix(
            self.dlist, self.plist)

    def testConstruct(self):
        self.assertEqual(self.m.data, self.data)
        self.assertEqual(self.m.positions, self.pos)
        self.assertEqual(self.m.nsam, 4)
        self.assertEqual(self.m.nsites, 2)

    def testModifyCppDataViaNumpy(self):
        d = np.array(self.m.data, copy=False)
        d[2] = 4
        x = [int(i) for i in self.m.site(0)]
        self.assertEqual(x[2], 4)


class testColumnViews(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.data = [0, 1, 1, 0, 0, 0, 0, 1]
        self.pos = [0.1, 0.2]
        self.m = libsequence.variant_matrix.VariantMatrix(self.data, self.pos)

    def testIterateColumns(self):
        for i in range(self.m.nsites):
            c = self.m.site(i)
            self.assertEqual(len(c), self.m.nsam)
        s = [j for j in self.m.site(0)]
        self.assertEqual(s, self.data[:self.m.nsam])
        s = [j for j in self.m.site(1)]
        self.assertEqual(s, self.data[self.m.nsam:])

    def testIterateRows(self):
        d = np.array(self.data, dtype=np.int8)
        d = d.reshape((self.m.nsites, self.m.nsam))
        for i in range(self.m.nsam):
            c = self.m.sample(i)
            s = np.array([j for j in c], dtype=np.int8)
            self.assertTrue(np.array_equal(s, d[:, i]))


if __name__ == "__main__":
    unittest.main()
