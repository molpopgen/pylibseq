import unittest
import pickle
from libsequence.polytable import *

class test_polytable(unittest.TestCase):
    def testNoPolyTable(self):
        with self.assertRaises(TypeError):
            x = PolyTable()
            
class test_simdata(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        d = [(0.1,"01010101"),(0.2,"01111111")]
        self.x = SimData(d)
    def testSimpleInit1(self):
        self.assertEqual(self.x.numsites(),2)
    def testSimpleInit2(self):
        self.assertEqual(self.x[0],b'00')
        self.assertEqual(self.x.size(),8)
    def testSimpleInit3(self):
        pos = [0.1,0.2]
        data = ["01","10"]
        x = SimData()
        x.assign(pos,data)
        self.assertEqual(x.numsites(),2)
    def testSimpleInit2(self):
        pos = [0.1,0.2]
        data = ["01","10"]
        x = SimData()
        x.assign(pos,data)
        self.assertEqual(x.size(),2)
    def testAssignFail1(self):
        with self.assertRaises(RuntimeError):
            ##Sample size at each site unequal
            d = [(0.1,"01010101"),(0.2,"0111")]
            x = SimData()
            x.assign(d)
    def testAssignFail2(self):
        with self.assertRaises(RuntimeError):
            pos = [0.1,0.2]
            #oops--3 sites in second haplotype
            data = ["01","100"]
            x = SimData()
            x.assign(pos,data)
    def testPickle(self):
        d = pickle.dumps(self.x,-1)
        x = pickle.loads(d)
        self.assertEqual(x.pos(),self.x.pos())
        self.assertEqual(x.data(),self.x.data())

class test_functions_PolySites(unittest.TestCase):
    def testRemoveMono1(self):
        d = [(0.1,"AGAG"),(0.2,"AAAA")]
        x = PolySites()
        x.assign(d)
        y=removeMono(x)
        self.assertEqual(y.numsites(),1)
        
class test_functions_SimData(unittest.TestCase):
    def testRemoveMono1(self):
        d = [(0.1,"01010101"),(0.2,"11111111")]
        x = SimData()
        x.assign(d)
        y=removeMono(x)
        self.assertEqual(y.numsites(),1)
    def testRemoveMono2(self):
        d = [(0.1,"00000000"),(0.2,"11111111")]
        x = SimData()
        x.assign(d)
        y=removeMono(x)
        self.assertEqual(y.numsites(),0)
    def testRemoveColumns1(self):
        d = [(0.1,"01010101"),
             (0.2,"11111111"),
             (0.3,"00010000"),
             (0.4,"00000001")
             ]
        x = SimData()
        x.assign(d)
        #Remove singletons with a lambda:
        y = removeColumns(x,lambda x:x.one > 1)
        pos = y.pos()
        self.assertEqual(pos,[0.1,0.2])
    def testIsValid(self):
        d = [(0.1,"01010101"),(0.2,"11011111")]
        x = SimData()
        x.assign(d)
        self.assertEqual(isValid(x),True)
    def testRemoveAmbig(self):
        ##Remove sites other than a,g,c,t,n,0,1
        d = [(0.1,"01010101"),(0.2,"1101Q111")]
        x = SimData()
        x.assign(d)
        y=removeAmbiguous(x)
        pos = y.pos()
        self.assertEqual(pos,[0.1])
    def testRemoveGaps(self):
        ##Remove sites other than a,g,c,t,n,0,1
        d = [(0.1,"01010101"),(0.2,"1101-111")]
        x = SimData()
        x.assign(d)
        y=removeGaps(x,gapchar='-')
        pos = y.pos()
        self.assertEqual(pos,[0.1])
                
if __name__ == '__main__':
    unittest.main()
