use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::str::FromStr;

use bio::io::fastq;
use bio::io::fastq::{Reader, Record, Records};

use crate::read_strategies::ten_x::TenXLayout;
use std::ops::Deref;
use std::borrow::Borrow;
use std::path::Path;

pub enum LayoutType {
    TENXV3,
    TENXV2,
    PAIREDUMI,
    PAIRED,
    SCI,
}

impl LayoutType {
    pub fn has_umi(&self) -> bool {
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

pub struct ReadLayout {
    read_one: Record,
    read_two: Option<Record>,
    index_one: Option<Record>,
    index_two: Option<Record>,
}

pub struct ReadIterator {
    pub read_one: Records<BufReader<File>>,
    pub read_two: Option<Records<BufReader<File>>>,
    pub index_one: Option<Records<BufReader<File>>>,
    pub index_two: Option<Records<BufReader<File>>>,
}

impl Iterator for ReadIterator {
    type Item = ReadLayout;

    fn next(&mut self) -> Option<ReadLayout> {
        let next_read_one = self.read_one.next();
        if next_read_one.is_some() {
            Some(ReadLayout {
                read_one: next_read_one.unwrap().unwrap(),
                read_two: match self.read_two
                {
                    Some(ref mut read_pointer) => Some(read_pointer.next().unwrap().unwrap()),
                    None => None,
                },
                index_one: match self.index_one
                {
                    Some(ref mut read_pointer) => Some(read_pointer.next().unwrap().unwrap()),
                    None => None,
                },
                index_two: match self.index_two
                {
                    Some(ref mut read_pointer) => Some(read_pointer.next().unwrap().unwrap()),
                    None => None,
                },
            })
        } else {
            None
        }
    }
}

impl ReadIterator
{
    pub fn new(read_1: String,
               read_2: String,
               index_1: String,
               index_2: String,
    ) -> ReadIterator  {
        ReadIterator {
            read_one: ReadIterator::open_reader(read_1).unwrap(),
            read_two: ReadIterator::open_reader(read_2),
            index_one: ReadIterator::open_reader(index_1),
            index_two: ReadIterator::open_reader(index_2),
        }
    }

    fn open_reader(filename: String) -> Option<Records<BufReader<File>>> {
        let check_path = Path::new(&filename).exists();
        if check_path {
            let f2 = File::open(filename.clone());
            let mut f2gz = fastq::Reader::new(File::open(filename).unwrap());

            let records = f2gz.records();
            Some(records)
        } else {
            None
        }
    }
}


pub trait SequenceLayout {
    fn name(&self) -> &Vec<u8>;
    fn umi(&self) -> Option<&Vec<u8>>;
    fn static_id(&self) -> Option<&Vec<u8>>;
    fn read_one(&self) -> &Vec<u8>;
    fn read_two(&self) -> Option<&Vec<u8>>;
    fn cell_id(&self) -> Option<&Vec<u8>>;
    fn layout_type(&self) -> LayoutType;
    fn get_unique_sequences(&self) -> Option<Vec<&Vec<u8>>>;
}


pub fn transform(read: ReadLayout, layout: &LayoutType) -> Box<dyn SequenceLayout> {
    match layout {
        LayoutType::TENXV3 => {
            assert!(read.read_two.is_some(), "Read two (read ID and UMI) must be defined for 10X");
            assert!(!read.index_one.is_some(), "Index 1 is invalid for 10X data");
            assert!(read.index_two.is_some(), "Index 2 is invalid for 10X data");

            let read2 = read.read_two.unwrap();
            let cell_id_sliced = read2.seq()[0..16].to_vec();
            let umi_sliced = read2.seq()[16..28].to_vec();

            Box::new(TenXLayout {
                name: read.read_one.id().as_bytes().to_vec(),
                umi: umi_sliced,
                cell_id: cell_id_sliced,
                read_one: read.read_one.seq().to_vec(),
            })
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