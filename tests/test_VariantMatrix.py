import unittest
import libsequence.variant_matrix


class testVariantMatrix(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.data = [0, 1, 1, 0, 0, 0, 0, 1]
        self.pos = [0.1, 0.2]

    def testConstruct(self):
        m = libsequence.variant_matrix.VariantMatrix(self.data, self.pos)
        self.assertEqual(m.data, self.data)
        self.assertEqual(m.positions, self.pos)
        self.assertEqual(m.nsam, 4)
        self.assertEqual(m.nsites, 2)


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
        import numpy as np
        d = np.array(self.data, dtype=np.int8)
        d = d.reshape((self.m.nsites, self.m.nsam))
        for i in range(self.m.nsam):
            c = self.m.sample(i)
            s = np.array([j for j in c], dtype=np.int8)
            self.assertTrue(np.array_equal(s, d[:, i]))


if __name__ == "__main__":
    unittest.main()
