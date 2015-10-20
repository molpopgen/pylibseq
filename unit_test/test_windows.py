import unittest

from libsequence.polytable import *
from libsequence.windows import *

class test_simData(unittest.TestCase):
    def testNumWindows(self):
        d = [(0.05,"01010101"),
             (0.1,"01010101"),
             (0.15,"01010101"),
             (0.2,"01010101"),
             (0.225,"01010101"),
             (0.25,"01010101"),
             (0.5,"01010101"),
             (0.95,"01010101"),
             (1.0,"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = simDataWindows(x,0.1,0.05,0,1)
        self.assertEqual(len(w),20)
    def testsWindowType(self):
        d = [(0.05,"01010101"),
             (0.1,"01010101"),
             (0.15,"01010101"),
             (0.2,"01010101"),
             (0.225,"01010101"),
             (0.25,"01010101"),
             (0.5,"01010101"),
             (0.95,"01010101"),
             (1.0,"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = simDataWindows(x,0.1,0.05,0,1)
        self.assertEqual(isinstance(w[0],simData),True)
    def testWindowContents(self):
        d = [(0.05,"01010101"),
             (0.1,"01010101"),
             (0.15,"01010101"),
             (0.2,"01010101"),
             (0.225,"01010101"),
             (0.25,"01010101"),
             (0.5,"01010101"),
             (0.95,"01010101"),
             (1.0,"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = simDataWindows(x,0.1,0.05,0,1)
        for i in range(len(w)):
            ##These are the positions that SHOULD be in each window...
            pp = [j[0] for j in d if (j[0] >= float(i)*0.05 and j[0] <= float(i)*0.05+0.1)]
            ##...but are they?
            ppwi = w[i].pos()
            self.assertEqual(ppwi,pp)
    def testWindowContents2(self):
        ##Same idea as previous test, but wacko window size/step len
        d = [(0.05,"01010101"),
             (0.1,"01010101"),
             (0.15,"01010101"),
             (0.2,"01010101"),
             (0.225,"01010101"),
             (0.25,"01010101"),
             (0.5,"01010101"),
             (0.95,"01010101"),
             (1.0,"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        stepsize = 0.032512
        wsize = 0.1312512
        w = simDataWindows(x,wsize,stepsize,0,1)
        for i in range(len(w)):
            ##These are the positions that SHOULD be in each window...
            pp = [j[0] for j in d if (j[0] >= float(i)*stepsize and j[0] <= float(i)*stepsize+wsize)]
            ##...but are they?
            ppwi = w[i].pos()
            self.assertEqual(ppwi,pp)
    def testException1(self):
        ##Fail with window size of 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,"01010101"),
                (0.1,"01010101"),
                (0.15,"01010101"),
                (0.2,"01010101"),
                (0.225,"01010101"),
                (0.25,"01010101"),
                (0.5,"01010101"),
                (0.95,"01010101"),
                (1.0,"01010101")]
            x = simData()
            x.assign(d)
            w = simDataWindows(x,0.,0.05,0,1)
    def testException2(self):
        ##Fail with step size of 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,"01010101"),
                (0.1,"01010101"),
                (0.15,"01010101"),
                (0.2,"01010101"),
                (0.225,"01010101"),
                (0.25,"01010101"),
                (0.5,"01010101"),
                (0.95,"01010101"),
                (1.0,"01010101")]
            x = simData()
            x.assign(d)
            w = simDataWindows(x,0.1,0,0,1)
    def testException3(self):
        ##Fail with window size < 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,"01010101"),
                (0.1,"01010101"),
                (0.15,"01010101"),
                (0.2,"01010101"),
                (0.225,"01010101"),
                (0.25,"01010101"),
                (0.5,"01010101"),
                (0.95,"01010101"),
                (1.0,"01010101")]
            x = simData()
            x.assign(d)
            w = simDataWindows(x,-1,0.05,0,1)
    def testException4(self):
        ##Fail with step size < 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,"01010101"),
                (0.1,"01010101"),
                (0.15,"01010101"),
                (0.2,"01010101"),
                (0.225,"01010101"),
                (0.25,"01010101"),
                (0.5,"01010101"),
                (0.95,"01010101"),
                (1.0,"01010101")]
            x = simData()
            x.assign(d)
            w = simDataWindows(x,0.1,-1,0,1)
if __name__ == '__main__':
    unittest.main()
