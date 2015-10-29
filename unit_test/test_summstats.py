import unittest
from libsequence.polytable import *
from libsequence.summstats import *

class test_nSL(unittest.TestCase):
    #API check
    def test_basic(self):
        d = [(0.1,"01010101"),(0.2,"01111111"),
             (0.3,"00011101"),(0.4,"11111000"),
             (0.5,"01010101"),(0.6,"00001111")
             ]
        x = simData()
        x.assign(d)
        stats = nSLiHS(x)
        self.assertEqual(len(stats),x.numsites())
        
    
