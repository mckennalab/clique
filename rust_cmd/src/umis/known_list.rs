use std::collections::{HashMap};
use std::collections::btree_set::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

use log::{info};

use crate::umis::sequence_clustering::BestHits;

use crate::alignment::fasta_bit_encoding::{FastaBase, reverse_complement};
use crate::read_strategies::sequence_layout::{UMIConfiguration, UMIPadding};

pub struct KnownList {
    pub name: UMIConfiguration,
    pub known_list_map: HashMap<Vec<FastaBase>, BestHits>,
    pub known_list_subset: HashMap<Vec<FastaBase>, Vec<Vec<FastaBase>>>,
    pub known_list_subset_key_size: usize,
}

impl KnownList {

    pub fn read_known_list_file(umi_type: &UMIConfiguration, starting_nmer_size: &usize) -> KnownList {

        let filename = umi_type.file.clone().unwrap();
        let filename = filename.as_str();

        info!("Setting up known list reader for {}", filename);
        let min_max = KnownList::get_min_max_length(filename);

        let rev_comp = umi_type.reverse_complement_sequences.clone().unwrap_or(false);
        let btree = KnownList::create_b_tree(filename, &umi_type.pad, &rev_comp, &min_max.1);

        // now that the barcodes are read and in-order, we'll make a prefix lookup. This should speed up matching later on
        let mut prefix: Option<Vec<FastaBase>> = None;
        let mut container: Vec<Vec<FastaBase>> = Vec::new();

        let mut existing_mapping = HashMap::new();
        let mut known_list_subset: HashMap<Vec<FastaBase>, Vec<Vec<FastaBase>>> = HashMap::new();

        for bytes in &btree {
            if !prefix.is_some() { prefix = Some(bytes[0..*starting_nmer_size].to_vec()); }

            existing_mapping.insert(bytes.clone(), BestHits { hits: vec![bytes.clone()], distance: 0 });

            let first_x = bytes.clone()[0..*starting_nmer_size].to_vec();
            if FastaBase::edit_distance(&first_x, prefix.as_ref().unwrap()) > 0 {
                known_list_subset.insert(prefix.unwrap(), container.clone());
                container.clear();
                prefix = Some(first_x);
            }
            container.push(bytes.clone());
        }
        known_list_subset.insert(prefix.unwrap(), container.clone());

        KnownList {
            name: umi_type.clone(),
            known_list_map: existing_mapping,
            known_list_subset,
            known_list_subset_key_size: starting_nmer_size.clone(),
        }
    }

    pub fn create_b_tree(filename: &str, pad_dir: &Option<UMIPadding>, reverse_comp: &bool, pad_length: &usize) -> BTreeSet<Vec<FastaBase>> {

        let raw_reader = BufReader::new(File::open(filename).unwrap());
        let mut btree = BTreeSet::new();
        for line in raw_reader.lines() {
            let mut bytes = line.unwrap();
            if *reverse_comp {
                bytes = FastaBase::string(&reverse_complement(&FastaBase::from_string(&bytes)));
            }
            match pad_dir {
                None => {
                    btree.insert(FastaBase::from_string(&bytes));
                }
                Some(x) => {
                    if bytes.len() < *pad_length {
                        match x {
                            UMIPadding::Left => {
                                let mut padded_bytes = bytes.clone();
                                padded_bytes.insert_str(0, &"-".repeat(pad_length - bytes.len()));
                                btree.insert(FastaBase::from_string(&padded_bytes));
                            }
                            UMIPadding::Right => {
                                let mut padded_bytes = bytes.clone();
                                padded_bytes.push_str(&"-".repeat(pad_length - bytes.len()));
                                btree.insert(FastaBase::from_string(&padded_bytes));
                            }
                        }
                    } else if bytes.len() == *pad_length {
                        btree.insert(FastaBase::from_string(&bytes));
                    } else {
                        panic!("Barcode {} is longer than the specified padding length of {}", &bytes, *pad_length);
                    }
                }
            }
        }
        btree
    }

    pub fn get_min_max_length(filename: &str) -> (usize,usize) {
        let raw_reader = BufReader::new(File::open(filename).unwrap());
        let mut min = usize::MAX;
        let mut max = usize::MIN;
        for line in raw_reader.lines() {
            let bytes = line.unwrap();
            let len = bytes.len();
            if len < min { min = len; }
            if len > max { max = len; }
        }
        (min,max)
    }
}