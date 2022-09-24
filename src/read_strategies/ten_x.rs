use super::sequence_layout::*;
use crate::read_strategies::sequence_file_containers::ReadSetContainer;
use crate::sorters::sorter::SortStructure;
use std::path::{Path, PathBuf};
use bio::io::fastq::Record;

pub struct TenXLayout {
    name: Vec<u8>,
    myumi: Vec<u8>,
    cellid: Vec<u8>,
    read_one: Vec<u8>,
    original_reads: Option<ReadSetContainer>,
}

impl TenXLayout {
    pub fn new(read: ReadSetContainer) -> TenXLayout {
        assert!(read.read_two.is_some(), "Read two (read ID and UMI) must be defined for 10X");
        assert!(!read.index_one.is_some(), "Index 1 is invalid for 10X data");
        assert!(!read.index_two.is_some(), "Index 2 is invalid for 10X data");

        let read2 = read.read_two.as_ref().unwrap().clone();
        let cell_id_sliced = read2.seq()[0..16].to_vec();
        let umi_sliced = read2.seq()[16..28].to_vec();

        TenXLayout {
            name: read.read_one.id().as_bytes().to_vec(),
            myumi: umi_sliced,
            cellid: cell_id_sliced,
            read_one: read.read_one.seq().to_vec(),
            original_reads: Some(read.clone()),
        }
    }
}

impl SequenceLayout for TenXLayout {
    fn name(&self) -> &Vec<u8> {
        &self.name
    }

    fn umi(&self) -> Option<&Vec<u8>> {
        Some(&self.myumi)
    }

    fn static_id(&self) -> Option<&Vec<u8>> {
        None
    }

    fn read_one(&self) -> &Vec<u8> {
        &self.read_one
    }

    fn read_two(&self) -> Option<&Vec<u8>> {
        None
    }

    fn cell_id(&self) -> Option<&Vec<u8>> {
        Some(&self.cellid)
    }

    fn layout_type(&self) -> LayoutType {
        LayoutType::TENXV3
    }

    fn original_reads(&self) -> Option<ReadSetContainer> {
        self.original_reads.clone()
    }

    fn has_original_reads(&self) -> bool {
        self.original_reads.is_some()
    }

    fn correct_known_sequence(&mut self, new: &Vec<u8>) {
        assert_eq!(self.cellid.len(),new.len());
        self.cellid = new.clone();
        let mut new_read2_seq = new.clone();
        let old_read_2 = self.original_reads.as_ref().unwrap().read_two.as_ref().unwrap();
        new_read2_seq.append(&mut old_read_2.seq().clone()[16..28].to_vec());
        let mut new_read2 = Record::with_attrs(old_read_2.id().clone(), None, new_read2_seq.as_slice(), old_read_2.qual().clone());
        self.original_reads = Some(ReadSetContainer::new_from_read2(new_read2, &self.original_reads.as_ref().unwrap()));
    }
}
