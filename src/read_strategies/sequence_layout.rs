use std::borrow::Borrow;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use bio::io::fastq;
use bio::io::fastq::{Reader, Record, Records};

use crate::read_strategies::ten_x::TenXLayout;
use crate::sorters::sorter::ReadSortingOnDiskContainer;
use crate::read_strategies::sequence_structures::ReadSetContainer;
use crate::sorters::sorter::SortStructure;
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum LayoutType {
    TENXV3,
    TENXV2,
    PAIREDUMI,
    PAIRED,
    SCI,
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
    fn umi(&self) -> Option<&Vec<u8>>;
    fn static_id(&self) -> Option<&Vec<u8>>;
    fn read_one(&self) -> &Vec<u8>;
    fn read_two(&self) -> Option<&Vec<u8>>;
    fn cell_id(&self) -> Option<&Vec<u8>>;
    fn layout_type(&self) -> LayoutType;
    //fn umi_sorting_path(&self) -> Option<Vec<SortStructure>>; // let ih = IteratorHolder { iter: 0..10 };
    fn original_reads(&self) -> Option<ReadSetContainer>;
    fn has_original_reads(&self) -> bool;
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
            unimplemented!()
        }
        LayoutType::TENXV2 => {
            unimplemented!()
        }

    }
}