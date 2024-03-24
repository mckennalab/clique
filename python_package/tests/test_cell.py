import unittest
from clique import tenX
from clique.cell import CellManager, construct_cell_tag


class TenXSingleCellStatsTest(unittest.TestCase):

    def testLoad(self):

        tenx = tenX.TenXSingleCellStats("test_data/maryam_test_data/outs/","test_data/3M-february-2018.txt.gz", False)
        tenx_cell_id_list = tenx.get_passing_cell_ids(True)
        print(len(tenx_cell_id_list))

        #manager = CellManager("test_data/maryam_test_data/maryam_10x_2023_10_17_post_butcher.bam",["e0"],["e1","e2"])
        manager = CellManager("test_data/maryam_test_data/maryam_10x_first_10K_reads.bam", ["e0"], ["e1", "e2"])


        manager.add_known_cell_id_list(tenx_cell_id_list)

        has_transcriptome,no_overlap = manager.intersection()
        #self.assertEqual(has_transcriptome, 16342)
