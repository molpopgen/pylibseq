import unittest
from libsequence.polytable import *
from libsequence.summstats import *

class test_nSL(unittest.TestCase):
    #API check
    def test_basic(self):
        d = [(0.1,b"01010101"),(0.2,b"01111111"),
             (0.3,b"00011101"),(0.4,b"11111000"),
             (0.5,b"01010101"),(0.6,b"00001111")
             ]
        x = SimData()
        x.assign(d)
        stats = nSLiHS(x)
        self.assertEqual(len(stats),x.numsites())

if __name__ == '__main__':
    unittest.main()
        
    
