use std::fs::File;
use std::io;
use std::io::BufReader;
use std::str::FromStr;

use bio::io::fastq;
use bio::io::fastq::{Reader, Record, Records};

use crate::read_strategies::ten_x::TenXLayout;
use std::ops::Deref;
use std::borrow::Borrow;

pub enum LayoutType {
    TENX_V3,
    TENX_V2,
    PAIREDUMI,
    PAIRED,
    SCI,
}

impl FromStr for LayoutType {
    type Err = ();

    fn from_str(input: &str) -> Result<LayoutType, Self::Err> {
        match input.to_uppercase().as_str() {
            "TENX_V3" => Ok(LayoutType::TENX_V3),
            "TENX_V2" => Ok(LayoutType::TENX_V2),
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
               read_2: Option<String>,
               index_1: Option<String>,
               index_2: Option<String>,
    ) -> ReadIterator  {
        ReadIterator {
            read_one: ReadIterator::open_reader(read_1).unwrap(),
            read_two: if read_2.is_some() {
                Some(ReadIterator::open_reader(read_2.unwrap()).unwrap())
            } else {
                None
            },
            index_one: if index_1.is_some() {
                Some(ReadIterator::open_reader(index_1.unwrap()).expect("Unable to open index one file").into_iter())
            } else {
                None
            },
            index_two: if index_2.is_some() {
                Some(ReadIterator::open_reader(index_2.unwrap()).expect("Unable to open index two file").into_iter())
            } else {
                None
            },
        }
    }

    fn open_reader(filename: String) -> Option<Records<BufReader<File>>> {
        let f2 = File::open(filename.clone());
        if f2.is_ok() {
            let mut f2gz = fastq::Reader::new(File::open(filename).unwrap());

            let records = f2gz.records();
            Some(records)
        } else {
            None
        }
    }
}


pub(crate) trait SequenceLayout {
    fn umi(&self) -> Option<&Vec<u8>>;
    fn static_id(&self) -> Option<&Vec<u8>>;
    fn read_one(&self) -> &Vec<u8>;
    fn read_two(&self) -> Option<&Vec<u8>>;
    fn cell_id(&self) -> Option<&Vec<u8>>;
    fn layout_type(&self) -> LayoutType;
}

fn transform(read_one: &Vec<u8>, read_two: Option<&Vec<u8>>, index_1: Option<&Vec<u8>>, index_2: Option<&Vec<u8>>, layout: LayoutType) -> Box<dyn SequenceLayout> {
    match layout {
        LayoutType::TENX_V3 => {
            assert!(read_two.is_some(), "Read two (read ID and UMI) must be defined for 10X");
            assert!(!index_1.is_some(), "Index 1 is invalid for 10X data");
            assert!(!index_2.is_some(), "Index 2 is invalid for 10X data");

            let cell_id_sliced = read_two.unwrap()[0..16].to_vec();
            let umi_sliced = read_two.unwrap()[16..28].to_vec();

            Box::new(TenXLayout {
                umi: umi_sliced,
                cell_id: cell_id_sliced,
                read_one: read_one.clone(),
            })
        }
        LayoutType::PAIREDUMI => {
            unimplemented!()
        }
        LayoutType::PAIRED => {
            unimplemented!()
        }
        LayoutType::SCI => {
            unimplemented!()
        }
        LayoutType::TENX_V2 => {
            unimplemented!()
        }
    }
}