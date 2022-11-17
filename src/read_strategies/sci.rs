use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::str;
use bio::io::fastq::Record;

use crate::read_strategies::sequence_file_containers::ReadSetContainer;
use crate::sorters::sorter::SortStructure;

use super::sequence_layout::*;

pub struct SciLayout {
    name: Vec<u8>,
    myumis: Box<HashMap<UMIType, Vec<u8>>>,
    read_one: Vec<u8>,
    read_two: Vec<u8>,
    original_reads: Option<ReadSetContainer>,
}

impl SciLayout {
    pub fn new(read: ReadSetContainer) -> SciLayout {
        assert!(read.read_two.is_some(), "Read two (read ID and UMI) must be defined for 10X");
        assert!(read.index_one.is_some(), "Index one (read ID and UMI) must be defined for 10X");
        assert!(!read.index_two.is_some(), "Index 2 is invalid for SCI data");


        let (ligation,umi_sliced,rt_index,static_id) = match (str::from_utf8(&read.read_one.seq()[9..14]).unwrap(), str::from_utf8(&read.read_one.seq()[10..15]).unwrap()) {
            (x, y) if y == "CTGTG" => {
                let ligation = read.read_one.seq()[0..10].to_vec();
                let umi = read.read_one.seq()[15..23].to_vec();
                let rt_index = read.read_one.seq()[23..33].to_vec();
                let static_id = read.read_one.seq()[58..68].to_vec();
                (ligation,umi,rt_index,static_id)
            }
            (_, _) => {
                let ligation = read.read_one.seq()[0..9].to_vec();
                let umi = read.read_one.seq()[14..22].to_vec();
                let rt_index = read.read_one.seq()[22..32].to_vec();
                let static_id = read.read_one.seq()[58..68].to_vec();
                (ligation,umi,rt_index,static_id)
            }
        };

        let pcr_index = read.index_two.clone().unwrap().clone();

        let remaining_read1 = read.read_one.seq()[68..read.read_one.seq().len()].to_vec();
        let read2 = read.read_two.as_ref().unwrap().clone().seq().to_vec();

        let umis: Box<HashMap<UMIType, Vec<u8>>> = Box::new(HashMap::from([(UMIType::SCILIG{size:ligation.len()}, ligation),
                                                                              (UMIType::SCIPCR{size:pcr_index.seq().len()}, pcr_index.seq().to_vec()),
                                                                              (UMIType::SCIRT{size:rt_index.len()}, rt_index),
            (UMIType::DEGENERATESEQ{size:static_id.len()}, static_id)]));

        SciLayout {
            name: read.read_one.id().as_bytes().to_vec(),
            myumis: umis,
            read_one: remaining_read1,
            read_two: read2,
            original_reads: Some(read.clone()),
        }
    }
}

impl SequenceLayout for SciLayout {
    fn name(&self) -> &Vec<u8> {
        &self.name
    }

    fn umis(&self) -> Option<&HashMap<UMIType, Vec<u8>>> {Some(&self.myumis)}

    fn read_one(&self) -> &Vec<u8> {
        &self.read_one
    }

    fn read_two(&self) -> Option<&Vec<u8>> {
        None
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

    fn correct_known_sequence(&mut self, new_type: UMIType, new_seq: &Vec<u8>) {
        match new_type {
            UMIType::SCIRT{size} => {
                self.myumis.insert(UMIType::SCIRT{size},new_seq.clone());
            },
            UMIType::SCILIG{size} => {
                self.myumis.insert(UMIType::SCILIG{size},new_seq.clone());
            },
            UMIType::SCIPCR{size} => {
                self.myumis.insert(UMIType::SCIPCR{size},new_seq.clone());
            },
            _ => panic!("Unsupported correction type {}",new_type),
        }
    }
}
