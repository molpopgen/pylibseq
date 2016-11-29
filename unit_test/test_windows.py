import unittest

from libsequence.polytable import *
from libsequence.windows import *

class test_simData(unittest.TestCase):
    def testNumWindows(self):
        d = [(0.05,b"01010101"),
             (0.1,b"01010101"),
             (0.15,b"01010101"),
             (0.2,b"01010101"),
             (0.225,b"01010101"),
             (0.25,b"01010101"),
             (0.5,b"01010101"),
             (0.95,b"01010101"),
             (1.0,b"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = Windows(x,0.1,0.05,0,1)
        self.assertEqual(len(w),20)
    def testsWindowType(self):
        d = [(0.05,b"01010101"),
             (0.1,b"01010101"),
             (0.15,b"01010101"),
             (0.2,b"01010101"),
             (0.225,b"01010101"),
             (0.25,b"01010101"),
             (0.5,b"01010101"),
             (0.95,b"01010101"),
             (1.0,b"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = Windows(x,0.1,0.05,0,1)
        self.assertEqual(isinstance(w[0],simData),True)
    def testWindowContents(self):
        d = [(0.05,b"01010101"),
             (0.1,b"01010101"),
             (0.15,b"01010101"),
             (0.2,b"01010101"),
             (0.225,b"01010101"),
             (0.25,b"01010101"),
             (0.5,b"01010101"),
             (0.95,b"01010101"),
             (1.0,b"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = Windows(x,0.1,0.05,0,1)
        for i in range(len(w)):
            ##These are the positions that SHOULD be in each window...
            pp = [j[0] for j in d if (j[0] >= float(i)*0.05 and j[0] <= float(i)*0.05+0.1)]
            ##...but are they?
            ppwi = w[i].pos()
            self.assertEqual(ppwi,pp)
    def testWindowContents2(self):
        ##Same idea as previous test, but wacko window size/step len
        d = [(0.05,b"01010101"),
             (0.1,b"01010101"),
             (0.15,b"01010101"),
             (0.2,b"01010101"),
             (0.225,b"01010101"),
             (0.25,b"01010101"),
             (0.5,b"01010101"),
             (0.95,b"01010101"),
             (1.0,b"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        stepsize = 0.032512
        wsize = 0.1312512
        w = Windows(x,wsize,stepsize,0,1)
        for i in range(len(w)):
            ##These are the positions that SHOULD be in each window...
            pp = [j[0] for j in d if (j[0] >= float(i)*stepsize and j[0] <= float(i)*stepsize+wsize)]
            ##...but are they?
            ppwi = w[i].pos()
            self.assertEqual(ppwi,pp)
    def testWindowContents3(self):
        ##Same as 1, but test iterability (?) of object
        d = [(0.05,b"01010101"),
             (0.1,b"01010101"),
             (0.15,b"01010101"),
             (0.2,b"01010101"),
             (0.225,b"01010101"),
             (0.25,b"01010101"),
             (0.5,b"01010101"),
             (0.95,b"01010101"),
             (1.0,b"01010101")]
        x = simData()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = Windows(x,0.1,0.05,0,1)
        i=0
        for wi in w:
            ##These are the positions that SHOULD be in each window...
            pp = [j[0] for j in d if (j[0] >= float(i)*0.05 and j[0] <= float(i)*0.05+0.1)]
            ##...but are they?
            ppwi = wi.pos()
            i+=1
            self.assertEqual(ppwi,pp)

    def testException1(self):
        ##Fail with window size of 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,b"01010101"),
                (0.1,b"01010101"),
                (0.15,b"01010101"),
                (0.2,b"01010101"),
                (0.225,b"01010101"),
                (0.25,b"01010101"),
                (0.5,b"01010101"),
                (0.95,b"01010101"),
                (1.0,b"01010101")]
            x = simData()
            x.assign(d)
            w = Windows(x,0.,0.05,0,1)
    def testException2(self):
        ##Fail with step size of 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,b"01010101"),
                (0.1,b"01010101"),
                (0.15,b"01010101"),
                (0.2,b"01010101"),
                (0.225,b"01010101"),
                (0.25,b"01010101"),
                (0.5,b"01010101"),
                (0.95,b"01010101"),
                (1.0,b"01010101")]
            x = simData()
            x.assign(d)
            w = Windows(x,0.1,0,0,1)
    def testException3(self):
        ##Fail with window size < 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,b"01010101"),
                (0.1,b"01010101"),
                (0.15,b"01010101"),
                (0.2,b"01010101"),
                (0.225,b"01010101"),
                (0.25,b"01010101"),
                (0.5,b"01010101"),
                (0.95,b"01010101"),
                (1.0,b"01010101")]
            x = simData()
            x.assign(d)
            w = Windows(x,-1,0.05,0,1)
    def testException4(self):
        ##Fail with step size < 0
        with self.assertRaises(RuntimeError):
            d = [(0.05,b"01010101"),
                (0.1,b"01010101"),
                (0.15,b"01010101"),
                (0.2,b"01010101"),
                (0.225,b"01010101"),
                (0.25,b"01010101"),
                (0.5,b"01010101"),
                (0.95,b"01010101"),
                (1.0,b"01010101")]
            x = simData()
            x.assign(d)
            w = Windows(x,0.1,-1,0,1)

class test_polyTable(unittest.TestCase):
    def testNumWindows(self):
        d = [(0.05,b"AGAGAGAG"),
             (0.1,b"AGAGAGAG"),
             (0.15,b"AGAGAGAG"),
             (0.2,b"AGAGAGAG"),
             (0.225,b"AGAGAGAG"),
             (0.25,b"AGAGAGAG"),
             (0.5,b"AGAGAGAG"),
             (0.95,b"AGAGAGAG"),
             (1.0,b"AGAGAGAG")]
        x = polySites()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = Windows(x,0.1,0.05,0,1)
        self.assertEqual(len(w),20)
    def testsWindowType(self):
        d = [(0.05,b"AGAGAGAG"),
             (0.1,b"AGAGAGAG"),
             (0.15,b"AGAGAGAG"),
             (0.2,b"AGAGAGAG"),
             (0.225,b"AGAGAGAG"),
             (0.25,b"AGAGAGAG"),
             (0.5,b"AGAGAGAG"),
             (0.95,b"AGAGAGAG"),
             (1.0,b"AGAGAGAG")]
        x = polySites()
        x.assign(d)
        ##We want sliding windows over the interval [0,1), step size = 0.1, jump size = 0.05
        ##There will be 20 such windows
        w = Windows(x,0.1,0.05,0,1)
        self.assertEqual(isinstance(w[0],polySites),True)
        
if __name__ == '__main__':
    unittest.main()
