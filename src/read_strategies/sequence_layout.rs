use std::borrow::Borrow;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::collections::HashMap;
use serde::{Serialize,Deserialize};
use bio::io::fastq;
use bio::io::fastq::{Reader, Record, Records};
use std::collections::BTreeMap;

use std::sync::{Arc, Mutex};


/// holds a set of reads for reading and writing to disk
#[derive(Serialize, Deserialize, Debug)]
pub struct ReadSetContainer {
    pub read_one: Record,
    pub read_two: Option<Record>,
    pub index_one: Option<Record>,
    pub index_two: Option<Record>,
}

impl Clone for ReadSetContainer {
    fn clone(&self) -> ReadSetContainer {
        ReadSetContainer {
            read_one: self.read_one.clone(),
            read_two: if self.read_two.as_ref().is_some() { Some(self.read_two.as_ref().unwrap().clone()) } else { None },
            index_one: if self.index_one.as_ref().is_some() { Some(self.index_one.as_ref().unwrap().clone()) } else { None },
            index_two: if self.index_two.as_ref().is_some() { Some(self.index_two.as_ref().unwrap().clone()) } else { None },
        }
    }
}

impl ReadSetContainer {
    pub fn new_from_read1(rec: Record) -> ReadSetContainer {
        ReadSetContainer {
            read_one: rec,
            read_two: None,
            index_one: None,
            index_two: None,
        }
    }
    pub fn new_from_read2(rec: Record, old: &ReadSetContainer) -> ReadSetContainer {
        ReadSetContainer {
            read_one: old.read_one.clone(),
            read_two: Some(rec),
            index_one: old.index_one.clone(),
            index_two: old.index_two.clone(),
        }
    }
}


#[derive(Clone)]
pub enum UMISortType {
    /// This enum represents how we should extract molecular identifiers and known sequences from the
    /// stream of reads. It specifies the location (as aligned to the reference), the maximum distance
    /// you allow before not considering two sequences to be from the same source.
    KNOWNTAG{name: String, ref_start: i32, ref_stop:i32, max_distance: usize},
    DEGENERATETAG{name: String, ref_start: i32, ref_stop:i32, max_distance: usize},
}


pub struct SequenceLayout {
    /// extracts read layout / sequences using known patterns taken from the YAML description file.
    name: Vec<u8>,
    umis: HashMap<UMISortType,Vec<u8>>,
    forward_seq: Vec<u8>,
    reverse_seq : Option<Vec<u8>>,
    original_reads: Option<ReadSetContainer>,
}

pub struct SequenceLayoutFactory {
    //sort_structure: Vec<SortStructure>,

}



impl SequenceLayoutFactory {
    pub fn from_yaml(yaml_file: String) {
        /// Load up a YAML document describing the layout of specific sequences within the reads. This configuration is specific
        /// to each sequencing platform and sequencing type (10X, sci, etc). The layout is described below.
        ///
        /// # Supported base tags
        ///
        /// *align*: (optional) * - contains one optional member, which can be _true_ or _false_
        /// *reads* - contains which read positions are required for this setup. The values, on individual lines, are _read1_, _read2_, _index1_, _index2_
        /// *umi* - contains one section per UMI, each with:
        /// - *name* - the base of each UMI section
        ///   - *read* - which read we can extract this from
        ///   - *start* - the starting position, if align = true is set this this is in relation to the reference, otherwise the offset into the read
        ///   - *length* - how long this sequence is
        ///   - *file* - (optional) which file contains known sequences that we should match to. One sequence per line, no header
        ///
        let mut file = File::open(&yaml_file).expect(&format!("Unable to open YAML configuration file: {}",&yaml_file));

        let mut yaml_contents = String::new();

        file.read_to_string(&mut yaml_contents)
            .expect(&format!("Unable to read contents of YAML configuration file: {}",&yaml_file));

        let docs = YamlLoader::load_from_str(&yaml_contents).unwrap();

        assert_eq!(docs.len(),1);

        let align: bool = match &docs[0]["align"].as_bool() {
            Some(x) => x.clone(),
            None => false
        };

        for doc_item in docs {

        }

    }
}
#[derive(Debug, PartialEq, Serialize, Deserialize)]
enum ReadPosition {
    READ1,
    READ2,
    INDEX1,
    INDEX2
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct UMIConfiguration {
    name: String,
    start: usize,
    length: usize,
    file: Option<String>,
}

impl UMIConfiguration {
    // pub fn from()
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct SequenceLayoutDesign {
    align: bool,
    reads: Vec<ReadPosition>,
    umi_configurations: Vec<UMIConfiguration>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str;
    use crate::utils::read_utils::fake_reads;


    #[test]
    fn test_basic_yaml_readback() {
        SequenceLayoutFactory::from_yaml(String::from("test_data/test_layout.yaml"))
    }
}
