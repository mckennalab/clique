import pysam


def construct_cell_tag(tags_to_values):
    class CellTags(object):
        __slots__ = tags_to_values.keys()
        def __init__(self, tags_to_values):
            for slot, arg in zip(CellTags.__slots__, tags_to_values.keys()):
                setattr(self, slot, tags_to_values[slot])
    return CellTags


class Cell:
    def __init__(self):

        self.barcode_sequences = []

    def add_barcodes(self, keys_and_values):
            self.barcode_sequences.append(construct_cell_tag(keys_and_values))

    def __repr__(self):
        return f"Cell with {len(self.barcode_sequences)} barcodes"

class CellManager:
    def __init__(self, bam_file_path, tags_that_define_cell, other_tags):
        self.cells = {}
        self.bam_file_path = bam_file_path
        self.tags_that_define_cell = tags_that_define_cell
        self.other_tags = other_tags
        self.process_bam_file()
        self.transcriptome_known_cell_ids = {}

    def add_known_cell_id_list(self,cell_id_list):
        for id in cell_id_list:
            self.transcriptome_known_cell_ids[id] = True

    def intersection(self):
        has_matching_transcriptome = 0
        missing_tanscriptome = 0
        for cell,tags in self.cells.items():
            if cell in self.transcriptome_known_cell_ids:
                has_matching_transcriptome += 1
            else:
                missing_tanscriptome += 1
        return(has_matching_transcriptome,missing_tanscriptome)

    def process_bam_file(self):

        with pysam.AlignmentFile(self.bam_file_path, "rb") as bam_file:
            for read in bam_file.fetch():
                tag_values = dict([(tag,read.get_tag(tag)) for tag in self.tags_that_define_cell])
                address = ".".join([x[1] for x in tag_values.items()])

                if address not in self.cells:
                    self.cells[address] = Cell()


                for tag in self.other_tags:
                    tag_values[tag] = read.get_tag(tag)

                self.cells[address].add_barcodes(tag_values)

    def get_cell(self, tag_values):
        return self.cells.get(tag_values, None)

    def __repr__(self):
        return f"CellManager with {len(self.cells)} cells"
