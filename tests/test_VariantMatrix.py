import unittest
import libsequence
import numpy as np


def is_singleton(x):
    from collections import Counter
    c = Counter([i for i in x])
    if len(c) == 2:
        for i in c.most_common():
            if i[1] == 1:
                return True
    return False


class testVariantMatrix(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.data = [0, 1, 1, 0, 0, 0, 0, 1]
        self.pos = [0.1, 0.2]
        self.m = libsequence.VariantMatrix(self.data, self.pos)

    def testConstruct(self):
        self.assertTrue(np.array_equal(self.m.data, np.array(
            self.data).reshape(self.m.nsites, self.m.nsam)))
        self.assertTrue(np.array_equal(self.m.positions, self.pos))
        self.assertEqual(self.m.nsam, 4)
        self.assertEqual(self.m.nsites, 2)

    def testModifyCppDataViaNumpy(self):
        d = np.array(self.m.data, copy=False)
        with self.assertRaises(ValueError):
            d[0][2] = 4

    def testFilterSites(self):
        self.assertEqual(self.m.nsites, 2)
        libsequence.filter_sites(
            self.m, is_singleton)
        self.assertEqual(self.m.nsites, 1)

    def testPickle(self):
        import pickle
        p = pickle.dumps(self.m, -1)
        up = pickle.loads(p)
        self.assertTrue(np.array_equal(up.data, self.m.data))
        self.assertTrue(np.array_equal(up.positions, self.m.positions))
        self.assertEqual(up.nsam, self.m.nsam)
        self.assertEqual(up.nsites, self.m.nsites)


class testColumnViews(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.data = [0, 1, 1, 0, 0, 0, 0, 1]
        self.pos = [0.1, 0.2]
        self.m = libsequence.VariantMatrix(self.data, self.pos)

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


class testCreationFromNumpy(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.d = np.array([0, 1, 1, 0, 0, 0, 0, 1],
                          dtype=np.int8).reshape((2, 4))
        self.pos = np.array([0.1, 0.2])

    def testConstruct(self):
        m = libsequence.VariantMatrix(self.d, self.pos)
        self.assertTrue(m.data.flags.writeable is False)
        # d is changed to non-writeable, too!
        self.assertTrue(self.d.flags.writeable is True)
        ma = np.array(m.data, copy=True)
        self.assertTrue(np.array_equal(
            np.sum(ma, axis=0), np.sum(self.d, axis=0)))
        self.assertTrue(np.array_equal(
            np.sum(ma, axis=1), np.sum(self.d, axis=1)))

    def testIteration(self):
        m = libsequence.VariantMatrix(self.d, self.pos)
        for i in range(m.nsites):
            c = m.site(i)
            self.assertEqual(len(c), m.nsam)
        s = np.array([j for j in m.site(0)])
        self.assertTrue(np.array_equal(s, self.d[0, ]))
        s = np.array([j for j in m.site(1)])
        self.assertTrue(np.array_equal(s, self.d[1, ]))

        for i in range(m.nsam):
            c = m.sample(i)
            s = np.array([j for j in c], dtype=np.int8)
            self.assertTrue(np.array_equal(s, self.d[:, i]))

    def testFilterSites(self):
        m = libsequence.VariantMatrix(self.d, self.pos)
        libsequence.filter_sites(m, is_singleton)


class testDataFromMsprime(unittest.TestCase):
    def testDirectConversionOfData(self):
        """
        Test is only run is msprime is
        available for import
        """
        try:
            import msprime
            try:
                ts = msprime.simulate(10, mutation_rate=10, random_seed=666)
                gm = ts.genotype_matrix()
                pos = np.array([i.position for i in ts.sites()])
                m = libsequence.VariantMatrix(gm, pos)
                # If the data conversion is correct, the row and
                # column sums must match
                ma = np.array(m, copy=False)
                self.assertTrue(np.array_equal(
                    np.sum(ma, axis=0), np.sum(gm, axis=0)))
                self.assertTrue(np.array_equal(
                    np.sum(ma, axis=1), np.sum(gm, axis=1)))
            except:
                raise
        except:
            pass


if __name__ == "__main__":
    unittest.main()
