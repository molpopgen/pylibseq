import unittest

from libsequence.polytable import simData
from libsequence.fst import fst

class test_fst(unittest.TestCase):
    def testException1(self):
        with self.assertRaises(RuntimeError):
            x = [(0.1,"0011"),(0.2,"1100"),
                (0.3,"0100"),(0.4,"1101"),
                (0.5,"0101")]
            d = simData()
            d.assign(x)
            ##the second argument's sum is > total sample size
            ##libsequence will throw a SeqException here,
            ##which gets tranlated to a RuntimeError
            f = fst(d,[2,3])
    def testShared(self):
        x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
        d = simData()
        d.assign(x)
        f = fst(d,[2,2])
        self.assertEqual(f.shared(0,1),[0.5])
    def testPriv(self):
        x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
        d = simData()
        d.assign(x)
        f = fst(d,[2,2])
        p = f.priv(0,1)
        expected = {0:[0.3],1:[0.4]}
        self.assertEqual(p,expected)
    def testFixed(self):
        x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
        d = simData()
        d.assign(x)
        f = fst(d,[2,2])
        p = f.fixed(0,1)
        expected = [0.1,0.2]
        self.assertEqual(p,expected)
        
if __name__ == '__main__':
    unittest.main()
