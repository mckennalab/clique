import os
import gzip
import unittest
from scipy.io import mmread


class TenXSingleCellStats:
    def __init__(self, ten_x_out_directory, matching_list, read_coverage):
        self.ten_x_out_directory = ten_x_out_directory
        print("Loading data from " + self.ten_x_out_directory)
        (self.filtered_list, self.unfiltered_list) = self.read10X_cell_lists()

        print("reading mapping file of barcodes")
        self.matching_list = {}
        self.map_feature_barcode(matching_list)

        if read_coverage:
            self.read_cell_coverage()

    def get_passing_cell_ids(self,mapped_to_known_tag):
        if mapped_to_known_tag:
            return [self.matching_list[x] for x in self.filtered_list]
        else:
            return self.filtered_list


    def map_feature_barcode(self, matching_list_file):
        '''
        Use the 10X matching capture ID list to map capture-tagged IDs to their corresponding cell IDs
        :param matching_list_file: the matching list provided by 10X, usually in lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
        :return:
        '''
        self.matching_list = {}
        match_file = gzip.open(matching_list_file, "rt")
        for line in match_file:
            tks = line.strip().split("\t")
            self.matching_list[tks[1]] = tks[0]

        self.filtered_list_matched = [self.matching_list[x] for x in self.filtered_list]

    def read10X_cell_lists(self):
        path_to_filtered_list = "filtered_feature_bc_matrix/barcodes.tsv.gz"
        unfiltered_list = "raw_feature_bc_matrix/barcodes.tsv.gz"

        filtered_list = read_10x_cell_list(os.path.join(self.ten_x_out_directory, path_to_filtered_list))
        unfiltered_list = read_10x_cell_list(os.path.join(self.ten_x_out_directory, unfiltered_list))

        # sanity check that the filtered list is fully a subset of the unfiltered list
        assert len(set(filtered_list).intersection(unfiltered_list)) == len(filtered_list)

        return ((filtered_list, unfiltered_list))

    def read_cell_coverage(self):
        assert hasattr(self, 'filtered_list')
        assert hasattr(self, 'path_to_unfiltered_list')
        matrix_market_file = "raw_feature_bc_matrix/matrix.mtx.gz"  # from the 'out' directory

        raw_coverage = mmread(os.path.join(self.ten_x_out_directory, matrix_market_file))
        self.unfiltered_cell_coverage = raw_coverage.sum(0)
        assert self.unfiltered_cell_coverage.shape[1] == len(self.unfiltered_list)

    def coverage_by_lineage_intersection(selfself, lineage_cell_ids):
        '''returns a data frame with the coverage of cells that do and don't have a corresponding lineage barcode'''
        assert hasattr(self, 'unfiltered_cell_coverage')


def read_10x_cell_list(cell_list_file):
    final_list = []
    print(cell_list_file)
    with gzip.open(cell_list_file, 'rt') as fl:
        for line in fl:
            final_list.append(line.split("-")[0])
    return final_list

