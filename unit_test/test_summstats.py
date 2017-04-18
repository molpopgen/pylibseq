import unittest
from libsequence.polytable import *
from libsequence.summstats import *

class test_nSL(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        d = [(0.1,"01010101"),(0.2,"01111111"),
             (0.3,"00011101"),(0.4,"11111000"),
             (0.5,"01010101"),(0.6,"00001111")
             ]
        self.x = SimData(d)
        
    #API check
    def test_nSLiHS(self):
        stats = nSLiHS(self.x)
        self.assertEqual(len(stats),self.x.numsites())
    def test_nSLiHS_gmap(self):
        gmap= {0.1:0.1,0.3:0.3,0.5:0.5,0.2:0.2,0.4:0.4,0.6:0.6}
        stats = nSLiHS(self.x,gmap)
        self.assertEqual(len(stats),self.x.numsites())

if __name__ == '__main__':
    unittest.main()
        
    
