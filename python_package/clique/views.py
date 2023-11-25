from collections.abc import Generator
from Levenshtein import distance
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
from clique import callers
import pysam
import numpy
import re


class CliqueReadSet:
    __slots__ = ["e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "name", "read", "alignment_rate",
                 "read_count", "cigar", "ref_diff_obj"]

    def __init__(self, tags, name, read, alignment_rate, ref_diff_obj, read_count, cigar):
        for tag_key, tag in tags.items():
            setattr(self, tag_key, tag)
        self.name = name
        self.read = read
        self.alignment_rate = alignment_rate
        self.read_count = read_count
        self.cigar = cigar
        self.ref_diff_obj = ref_diff_obj

    def reference_difference(self):
        return (self.ref_diff_obj.create_reference_difference_array(self.read, self.cigar))

    def __str__(self):
        return "-".join([str(getattr(self, x)) if hasattr(self, x) else '' for x in self.__slots__])


class LineageBamFile(Generator):
    def __init__(self, bam_file, reference, minimum_alignment_rate, minimum_read_count, anchors):
        self.bam_file_name = bam_file
        self.reads = []
        self.reference_raw = reference
        self.reference = EventCaller(ReferenceDifference.reference_number_masking(reference))
        self.anchors = anchors
        self.minimum_alignment_rate = minimum_alignment_rate
        self.minimum_read_count = minimum_read_count
        self.bamfile = pysam.AlignmentFile(self.bam_file_name, "rb").fetch(until_eof=True)
        self.reviewed_reads = 0

    def clone(self):
        return LineageBamFile(self.bam_file_name,
                              self.reference_raw,
                              self.minimum_alignment_rate,
                              self.minimum_read_count,
                              self.anchors)

    def send(self, ignored_arg):
        while True:
            self.reviewed_reads += 1
            read = self.bamfile.__next__()

            tags = {}
            alignment_rate = 0.0
            read_count = 0
            for tag in read.tags:
                if tag[0] >= 'e0' and tag[0] <= 'e9':
                    tags[tag[0]] = tag[1]
                if tag[0] == 'rm':
                    alignment_rate = float(tag[1])
                if tag[0] == 'rc':
                    read_count = int(tag[1])

            if alignment_rate >= self.minimum_alignment_rate and read_count >= self.minimum_read_count and sum(
                    [anchor in read.query_sequence for anchor in self.anchors]):
                return CliqueReadSet(tags, read.query_name, read.query_sequence, alignment_rate, self.reference,
                                     read_count, read.cigarstring)

    def throw(self, type=None, value=None, traceback=None):
        raise StopIteration


class BaseCalledCell:
    def __init__(self, cell_id):
        self.cell_id = cell_id
        self.integration_ids = []
        self.editing_outcomes = {}  # mapping integration IDs to outcomes
        self.read_counts = []

    def add_editing(self, integration_id, editing_outcomes, read_count):
        if integration_id in self.integration_ids:
            self.editing_outcomes[integration_id].append(editing_outcomes)
            self.read_counts[self.integration_ids.index(integration_id)] += read_count
        else:
            self.integration_ids.append(integration_id)
            self.editing_outcomes[integration_id] = [editing_outcomes]
            self.read_counts.append(read_count)



class CellList:

    def cluster_integration_ids(distances):

        condensed_dist_matrix = squareform(distances)

        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=2, metric="precomputed",
                                             linkage="single").fit(distances)
        clustering.labels_

        b = [int_list[index] for index, x in enumerate(clustering.labels_) if x == 9]
        id_to_clone = {int_list[index]: x for index, x in enumerate(clustering.labels_)}

    def ids_to_distances(list1, list2):
        overlap_count = 0
        distances = np.zeros((len(list1), len(list2)))
        for indexi, i in enumerate(list1):
            for indexj, j in enumerate(list2):
                distances[indexi, indexj] = distance(i, j)
                if distances[indexi, indexj] <= 1 and indexi < indexj:
                    print(i, " ", j)
                    overlap_count += 1
        print(overlap_count)
        print(len(list1))
        return distances
    def comp(e):
        if e == 'A':
            return 'T'
        elif e == 'C':
            return 'G'
        elif e == 'G':
            return 'C'
        elif e == 'T':
            return 'A'
        return 'N'

    def rev_comp(e):
        comp_str = [comp(x) for x in e]
        comp_str.reverse()
        return "".join(comp_str)

    #known_integration_list = []
    #known_integration_list_file = open(
    #    "/dartfs/rc/lab/M/McKennaLab/projects/Maryam/2022_11_30_MFN11/MFN11_data/clique/2023_09_27_full_barcode_mashup_filtered.txt")
    #known_integration_list_file.readline()
    #for line in known_integration_list_file:
    #    known_integration_list.append(rev_comp(line.split("\t")[1]))


    #int_list = check_known_integration_list_overlap(known_integration_list, msl)
    def check_known_integration_list_overlap(integration_list, becl):
        matching_cell_level = []
        integration_matching_level = {}
        all_integrations = {}

        for cell_id, cell in becl.matched_cells.items():

            matched_in_cell = 0
            for cid in cell.integration_ids:
                all_integrations[cid] = all_integrations.get(cid, 0) + 1

            for known_integration in integration_list:
                if known_integration in cell.integration_ids:
                    matched_in_cell += 1
                    integration_matching_level[known_integration] = integration_matching_level.get(known_integration,
                                                                                                   0) + 1
            matching_cell_level.append(matched_in_cell)

        return [k for k, v in all_integrations.items() if v > 4000]




class BaseEditingCellList:
    '''
    Collect the editing outcomes for a list of cells from a base-editing recorder sequencing library
    '''

    def __init__(self, bam_iterator, single_cell_obj, cell_id_tag, integration_id_tag):
        '''

        Keyword arguments:
        bam_iterator - a open bam file iterator used to parse out reads
        single_cell_obj - we use this object to match against the known cell IDs
        cell_id_tag - which tags we use to match single-cell transcripomes with
        integration_id_tag - the tag that identifies the integration ID within the cell
        '''
        self.single_cell_obj = single_cell_obj
        self.matched_cell_barcodes = 0
        self.unmatched_cell_barcodes = 0
        self.no_tag = 0
        self.matched_cells = {x: BaseCalledCell(x) for x in single_cell_obj.filtered_list_matched}

        for index, barcode_read in enumerate(bam_iterator):
            cell_id = getattr(barcode_read, cell_id_tag)
            if cell_id in self.matched_cells:
                self.matched_cells[cell_id].add_editing(getattr(barcode_read, integration_id_tag),
                                                        barcode_read.reference_difference(), barcode_read.read_count)
                self.matched_cell_barcodes += 1
            else:
                self.unmatched_cell_barcodes += 1

            if index % 10000000 == 0:
                print("Processed ", index, " reads")