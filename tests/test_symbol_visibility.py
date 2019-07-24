import unittest
import libsequence
import cpptests


class TestAlleleCountMatrix(unittest.TestCase):
    def test_create(self):
        am = cpptests.create_AlleleCountMatrix()
        self.assertTrue(isinstance(am, libsequence.AlleleCountMatrix))


if __name__ == "__main__":
    unittest.main()
