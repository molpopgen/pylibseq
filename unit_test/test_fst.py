import unittest

from libsequence.polytable import SimData
from libsequence.fst import Fst

class test_Fst(unittest.TestCase):
    def testException1(self):
        with self.assertRaises(RuntimeError):
            x = [(0.1,"0011"),(0.2,"1100"),
                (0.3,"0100"),(0.4,"1101"),
                (0.5,"0101")]
            d = SimData(x)
            ##the second argument's sum is > total sample size
            ##libsequence will throw a std::runtime_error here,
            ##which gets tranlated to a RuntimeError
            f = Fst(d,[2,3])
    def testShared(self):
        x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
        d = SimData(x)
        f = Fst(d,[2,2])
        self.assertEqual(f.shared(0,1),[0.5])
    def testPriv(self):
        x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
        d = SimData(x)
        f = Fst(d,[2,2])
        p = f.priv(0,1)
        expected = {0:[0.3],1:[0.4]}
        self.assertEqual(p,expected)
    def testFixed(self):
        x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
        d = SimData(x)
        f = Fst(d,[2,2])
        p = f.fixed(0,1)
        expected = [0.1,0.2]
        self.assertEqual(p,expected)
    #shared/fixed/private will get exceptions from libsequence
    #if i,j are out of range.
    def testExceptionShared(self):
        with self.assertRaises(ValueError):
            x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
            d = SimData(x)
            f = Fst(d,[2,2])
            #2 is out of range.
            sh = f.shared(2,1)
    def testExceptionFixed(self):
        with self.assertRaises(ValueError):
            x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
            d = SimData(x)
            f = Fst(d,[2,2])
            #2 is out of range.
            sh = f.fixed(2,1)
    def testExceptionPriv(self):
        with self.assertRaises(ValueError):
            x = [(0.1,"0011"),(0.2,"1100"),
            (0.3,"0100"),(0.4,"1101"),
            (0.5,"0101")]
            d = SimData(x)
            f = Fst(d,[2,2])
            #2 is out of range.
            sh = f.priv(2,1)
        
if __name__ == '__main__':
    unittest.main()
