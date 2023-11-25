#
#bam_name = "/home/Kiewit/f003w5r/McKennaLab/projects/clique/clique_versions/git_clique/rust_cmd/maryam_data/maryam_10x_2023_10_17_post_butcher.bam"
#lb_iterator = LineageBamFile(bam_name,
#                   "",
#                   0.90,
#                   1,
#                   "GTTGATAC")

import unittest
from clique import views
from clique import tenX
bam_file = "test_data/10000_base_edited_reads.bam"
reference = "0000000000000000111111111111aagcagtggtatcaacgcagagtacatgggCCAGGAAGTACTCGAGTACTTCCTGGacatgtCCTGTCATCTTAGCTAAGATGACAGGtacatc222222222222aacgtaannnnnnnnnnnnnnnnnnnnnnnnnnnnctgtctcttatacacatctgacgctgccgacgannnnnnnnnngtgtagatctcggtggtcgccgtatcattttacgtt"
matching_prop = 0.90
minimum_read_count = 1
anchors = "GTTGATAC"

#class LineageBamFileTest(unittest.TestCase):

    #def testLoad(self):

        #lineage_view = views.LineageBamFile(bam_file, reference, matching_prop, minimum_read_count, anchors)



    # def testBaseEditingCellList(self):

        #tenx = tenX.TenXSingleCellStats("test_data/outs/","test_data/3M-february-2018.txt.gz", False)

        #lineage_view = views.LineageBamFile(bam_file, reference, matching_prop, minimum_read_count, anchors)

        #base_editing = views.BaseEditingCellList(lineage_view, tenx, "e0", "e2")

    #if __name__ == '__main__':
    #    unittest.main()