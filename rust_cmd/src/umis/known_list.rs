use std::collections::{HashMap};
use std::collections::btree_set::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

use log::{info};

use crate::umis::sequence_clustering::BestHits;

use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::sequence_layout::UMIConfiguration;

pub struct KnownList {
    pub name: UMIConfiguration,
    pub known_list_map: HashMap<Vec<FastaBase>, BestHits>,
    pub known_list_subset: HashMap<Vec<FastaBase>, Vec<Vec<FastaBase>>>,
    pub known_list_subset_key_size: usize,
}

impl KnownList {

    pub fn read_known_list_file(umi_type: &UMIConfiguration, filename: &str, starting_nmer_size: &usize) -> KnownList {

        info!("Setting up known list reader for {}", &filename);
        let btree = KnownList::create_b_tree(filename);

        // now that the barcodes are read and in-order, we make a quick lookup prefix lookup to speed up matching later on
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

    pub fn create_b_tree(filename: &str) -> BTreeSet<Vec<FastaBase>> {

        let raw_reader = BufReader::new(File::open(filename).unwrap());


        // sort the barcodes as we go with a btree
        let mut btree = BTreeSet::new();

        for line in raw_reader.lines() {
            let bytes = line.unwrap();
            if bytes.contains("-") {
                btree.insert(FastaBase::from_vec_u8(&bytes.split("-").collect::<Vec<&str>>().get(0).unwrap().as_bytes().to_vec()));
            } else {
                btree.insert(FastaBase::from_string(&bytes));
            }
        }
        btree
    }
}