import unittest
from libsequence.polytable import *

class test_polytable(unittest.TestCase):
    def testNoPolyTable(self):
        with self.assertRaises(RuntimeError):
            x = polyTable()
            
class test_simdata(unittest.TestCase):
    def testSimpleInit1(self):
        d = [(0.1,"01010101"),(0.2,"01111111")]
        x = simData()
        x.assign(d)
        self.assertEqual(x.numsites(),2)
    def testSimpleInit2(self):
        d = [(0.1,"01010101"),(0.2,"01111111")]
        x = simData()
        x.assign(d)
        self.assertEqual(x.size(),8)
    def testSimpleInit3(self):
        pos = [0.1,0.2]
        data = ["01","10"]
        x = simData()
        x.assign_sep(pos,data)
        self.assertEqual(x.numsites(),2)
    def testSimpleInit2(self):
        pos = [0.1,0.2]
        data = ["01","10"]
        x = simData()
        x.assign_sep(pos,data)
        self.assertEqual(x.size(),2)
    def testAssignFail1(self):
        with self.assertRaises(RuntimeError):
            ##Sample size at each site unequal
            d = [(0.1,"01010101"),(0.2,"0111")]
            x = simData()
            x.assign(d)
    def testAssignFail2(self):
        with self.assertRaises(RuntimeError):
            pos = [0.1,0.2]
            #oops--3 sites in second haplotype
            data = ["01","100"]
            x = simData()
            x.assign_sep(pos,data)

class test_functions_simData(unittest.TestCase):
    def testRemoveMono1(self):
        d = [(0.1,"01010101"),(0.2,"11111111")]
        x = simData()
        x.assign(d)
        removeMono(x)
        self.assertEqual(x.numsites(),1)
    def testRemoveMono2(self):
        d = [(0.1,"00000000"),(0.2,"11111111")]
        x = simData()
        x.assign(d)
        removeMono(x)
        self.assertEqual(x.numsites(),0)
    def testFreqFilter1(self):
        d = [(0.1,"01010101"),
             (0.2,"11111111"),
             (0.3,"00010000"),
             (0.4,"00000001")
             ]
        x = simData()
        x.assign(d)
        #Remove singletons:
        freqFilter(x,2)
        pos = x.pos()
        self.assertEqual(pos,[0.1,0.2])
        
if __name__ == '__main__':
    unittest.main()
