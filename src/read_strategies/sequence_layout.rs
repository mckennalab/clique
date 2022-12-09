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
use crate::sorters::known_list::KnownList;
use std::sync::{Arc, Mutex};

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
    TENXRT,
    SCIRT,
    SCILIG,
    SCIPCR,
    DEGENERATESEQ
}

pub struct UMIInstance {
    pub umi_type: UMIType,
    pub umi_length: usize,
    pub umi_file: Option<String>,
    pub known_list: Option<Arc<Mutex<KnownList>>>,
}



impl FromStr for UMIType {
    type Err = ();

    fn from_str(input: &str) -> Result<UMIType, Self::Err> {

        match input {
            "TENXRT" => Ok(UMIType::TENXRT),
            "SCILIG" => Ok(UMIType::SCILIG),
            "SCIPCR" => Ok(UMIType::SCIPCR),
            "SCIRT" => Ok(UMIType::SCIRT),
            _ => Err(()),
        }
    }
}

impl std::fmt::Display for UMIType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            UMIType::TENXRT => {let res = write!(f, "TENXRT");res}
            UMIType::SCIRT => {let res = write!(f, "SCIRT");res}
            UMIType::SCILIG => {let res = write!(f, "SCILIG");res}
            UMIType::SCIPCR => {let res = write!(f, "SCIPCR");res}
            UMIType::DEGENERATESEQ => {let res = write!(f, "DEGENERATESEQ");res}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn umi_type_from_string() {
        let umi_string = String::from("SCIRT").as_bytes().to_owned();

        let type_of = UMIType::from_str(std::str::from_utf8(&umi_string).unwrap());

        assert!(type_of.is_ok(),"Unable to convert SCIRT into a UMIType");

    }

}