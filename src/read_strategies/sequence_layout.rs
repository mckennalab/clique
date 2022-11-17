use std::borrow::Borrow;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::collections::HashMap;

use bio::io::fastq;
use bio::io::fastq::{Reader, Record, Records};

use crate::read_strategies::ten_x::TenXLayout;
use crate::read_strategies::sequence_file_containers::ReadSetContainer;
use crate::sorters::sorter::SortStructure;

use crate::read_strategies::sequence_layout::LayoutType::SCI;
use crate::read_strategies::sci::SciLayout;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum LayoutType {
    TENXV3,
    TENXV2,
    PAIREDUMI,
    PAIRED,
    SCI,
}

impl std::fmt::Display for LayoutType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            LayoutType::TENXV3 => {let res = write!(f,"TENXV3");res}
            LayoutType::TENXV2 => {let res = write!(f,"TENXV2");res}
            LayoutType::PAIREDUMI => {let res = write!(f,"PAIREDUMI");res}
            LayoutType::PAIRED => {let res = write!(f,"PAIRED");res}
            LayoutType::SCI => {let res = write!(f,"SCI");res}
        }
    }
}

impl LayoutType {
    pub fn has_tags(&self) -> bool {
        match *self {
            LayoutType::TENXV3 => true,
            LayoutType::TENXV2 => true,
            LayoutType::PAIREDUMI => true,
            LayoutType::PAIRED => true,
            LayoutType::SCI => true,
        }
    }
}

impl FromStr for LayoutType {
    type Err = ();

    fn from_str(input: &str) -> Result<LayoutType, Self::Err> {
        match input.to_uppercase().as_str() {
            "TENXV3" => Ok(LayoutType::TENXV3),
            "TENXV2" => Ok(LayoutType::TENXV2),
            "PAIREDUMI" => Ok(LayoutType::PAIREDUMI),
            "PAIRED" => Ok(LayoutType::PAIRED),
            "SCI" => Ok(LayoutType::SCI),
            _ => Err(()),
        }
    }
}

pub struct SequenceLayoutFactory {
    pub known_list_file: HashMap<LayoutType,PathBuf>
}


pub trait SequenceLayout {
    fn name(&self) -> &Vec<u8>;
    fn umis(&self) -> Option<&HashMap<UMIType,Vec<u8>>>;
    fn read_one(&self) -> &Vec<u8>;
    fn read_two(&self) -> Option<&Vec<u8>>;
    fn layout_type(&self) -> LayoutType;
    fn original_reads(&self) -> Option<ReadSetContainer>;
    fn has_original_reads(&self) -> bool;
    fn correct_known_sequence(&mut self, umi_type: UMIType, new_seq: &Vec<u8>);
}

pub fn transform(read: ReadSetContainer, layout: &LayoutType) -> Box<dyn SequenceLayout> {
    match layout {
        LayoutType::TENXV3 => {
            Box::new(TenXLayout::new(read))
        }
        LayoutType::PAIREDUMI => {
            unimplemented!()
        }
        LayoutType::PAIRED=> {
            unimplemented!()
        }
        LayoutType::SCI => {
            Box::new(SciLayout::new(read))
        }
        LayoutType::TENXV2 => {
            unimplemented!()
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum UMIType {
    TENXRT{size: usize},
    SCIRT{size: usize},
    SCILIG{size: usize},
    SCIPCR{size: usize},
    DEGENERATESEQ{size: usize}
}


impl FromStr for UMIType {
    type Err = ();

    fn from_str(input: &str) -> Result<UMIType, Self::Err> {
        let tokens: Vec<&str> = input.split(' ').collect();
        if tokens.len() != 2 {
            return Err(())
        }
        match (tokens.get(0).unwrap().to_uppercase().as_str(),usize::from_str(&tokens.get(1).unwrap().to_uppercase()).unwrap()) {
            ("TENXRT",x) => Ok(UMIType::TENXRT{size: x}),
            ("SCIRT",x) => Ok(UMIType::SCIRT{size: x}),
            ("SCILIG",x) => Ok(UMIType::SCILIG{size: x}),
            ("SCIPCR",x) => Ok(UMIType::SCIPCR{size: x}),
            _ => Err(()),
        }
    }
}

impl std::fmt::Display for UMIType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            UMIType::TENXRT{size} => {let res = write!(f, "TENXRT:{}",size);res}
            UMIType::SCIRT{size} => {let res = write!(f, "SCIRT:{}",size);res}
            UMIType::SCILIG{size} => {let res = write!(f, "SCILIG:{}",size);res}
            UMIType::SCIPCR{size} => {let res = write!(f, "SCIPCR:{}",size);res}
            UMIType::DEGENERATESEQ { size } => {let res = write!(f, "DEGENERATESEQ:{}",size);res}
        }
    }
}