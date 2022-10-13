use std::collections::{HashMap};
use std::str;
use crate::knownlist::{KeyValueDualMapping, read_known_list_file};
use std::str::FromStr;
use crate::knownlist::read_lines;

pub enum BarcodeType {
    STATIC,
    DEGENERATE,
}

impl FromStr for BarcodeType {
    type Err = ();

    fn from_str(input: &str) -> Result<BarcodeType, Self::Err> {
        match input {
            "STATIC"  => Ok(BarcodeType::STATIC),
            "DEGENERATE"  => Ok(BarcodeType::DEGENERATE),
            _      => Err(()),
        }
    }
}


pub struct Configuration {
    pub symbol_to_codes: HashMap<String, KeyValueDualMapping>,
    pub symbol_to_mismatches: HashMap<String, u32>,
    pub symbol_to_name: HashMap<String, String>,
    pub symbol_to_type: HashMap<String, BarcodeType>,
}

pub fn read_configuration_file(filename: &String) -> Configuration {
    let mut symbol_to_codes = HashMap::new();
    let mut symbol_to_mismatches = HashMap::new();
    let mut symbol_to_name = HashMap::new();
    let mut symbol_to_type = HashMap::new();
    print!("Read file {}", filename);

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(filename) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                let split = ip.split("\t");
                let vec: Vec<&str> = split.collect();
                if vec.len() != 5 {
                    trace!("line: {:?} not in five parts (name,symbol,type,mismatches,file)", vec);
                }

                symbol_to_name.insert(vec[1].to_owned(),vec[0].to_owned());
                symbol_to_type.insert(vec[1].to_owned(),BarcodeType::from_str(vec[2]).unwrap());
                symbol_to_mismatches.insert(vec[1].to_owned(),u32::from_str_radix(vec[3],10).unwrap());

                let file = read_known_list_file(&vec[4].to_owned());
                if file.is_some() {
                    symbol_to_codes.insert(vec[1].to_owned(), file.unwrap());
                } else {
                    panic!("Unable to parse file {}",vec[4]);
                }
            }
        }
    }
    Configuration {
        symbol_to_codes,
        symbol_to_mismatches,
        symbol_to_name,
        symbol_to_type
    }
}