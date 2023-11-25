"""This module tests function reverse_words
provided by module mod.py."""
import unittest
from clique import tenX


class TenXSingleCellStatsTest(unittest.TestCase):

    def testLoad(self):

        tenx = tenX.TenXSingleCellStats("test_data/outs/","test_data/3M-february-2018.txt.gz", False)

        self.assertEqual(len(tenx.matching_list),6794880)

        self.assertEqual(len(tenx.filtered_list_matched),17742)

    if __name__ == '__main__':
        unittest.main()