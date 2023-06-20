use std::collections::{HashMap, VecDeque};
use std::collections::btree_set::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::slice::Iter;

use log::{info, warn};
use rayon::prelude::*;
use rayon::prelude::IntoParallelRefMutIterator;

use crate::fasta_comparisons::*;
use crate::umis::sequence_clustering::*;
use crate::umis::sequence_clustering::BestHits;
use crate::utils::base_utils::edit_distance;
use indicatif::ProgressBar;
use std::str::FromStr;
use std::num::ParseIntError;
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::sequence_layout::UMIConfiguration;


pub struct KnownListConsensus {
    observed_identifier_to_best_matches: HashMap<Vec<FastaBase>, BestHits>,
    observed_identifiers_counts: HashMap<Vec<FastaBase>, u64>,
    total_reads: u64,
}

impl KnownListConsensus {
    pub fn new() -> KnownListConsensus {
        let ids: HashMap<Vec<FastaBase>, BestHits> = HashMap::new();
        let counts: HashMap<Vec<FastaBase>, u64> = HashMap::new();
        KnownListConsensus { observed_identifier_to_best_matches: ids, observed_identifiers_counts: counts, total_reads: 0 }
    }

    pub fn add_hit(&mut self, barcode: &Vec<FastaBase>, best_hit: BestHits) {
        self.total_reads = self.total_reads + 1;
        let new_total = *self.observed_identifiers_counts.entry(barcode.clone()).or_insert(0) + 1 as u64;
        self.observed_identifiers_counts.insert(barcode.clone(), new_total);
        self.observed_identifier_to_best_matches.insert(barcode.clone(), best_hit);
    }


    pub fn match_to_knownlist(&self) -> KnownListConsensus {
        let mut total = 0;
        let mut unresolved = 0;
        let mut resolved = 0;
        let mut clean_mapping: HashMap<Vec<FastaBase>, BestHits> = HashMap::new();
        let mut clean_ids: HashMap<Vec<FastaBase>, u64> = HashMap::new();

        self.observed_identifier_to_best_matches.iter().for_each(|(k, v)| {
            let count = *self.observed_identifiers_counts.get(&k.clone()).unwrap();
            total += count;
            if v.hits.len() > 1 {
                let matching: Vec<&Vec<FastaBase>> = v.hits.iter().filter(|hit| self.observed_identifiers_counts.contains_key(*hit)).collect::<Vec<&Vec<FastaBase>>>();
                if matching.len() != 1 {
                    unresolved += count;
                } else {
                    resolved += count;
                    let best_hits = BestHits { hits: vec![matching[0].clone()], distance: v.distance };
                    clean_mapping.insert(k.clone(), best_hits);
                    clean_ids.insert(k.clone(), count);
                }
            } else {
                clean_mapping.insert(k.clone(), v.clone());
                clean_ids.insert(k.clone(), count);
            }
        });

        info!("total = {}, unresolved = {}, resolved = {}, final = {}", total, unresolved, resolved, clean_mapping.len());

        KnownListConsensus {
            observed_identifier_to_best_matches: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total,
        }
    }

    pub fn filter_to_minimum_coverage(&self, minimum_coverage: u64) -> KnownListConsensus {
        let mut clean_mapping: HashMap<Vec<FastaBase>, BestHits> = HashMap::new();
        let mut clean_ids: HashMap<Vec<FastaBase>, u64> = HashMap::new();

        let mut total_count = 0;

        self.observed_identifiers_counts.iter().map(|(k, count)| {
            if *count >= minimum_coverage {
                clean_ids.insert(k.clone(), *count);
                clean_mapping.insert(k.clone(), self.observed_identifier_to_best_matches.get(k).unwrap().clone());
                total_count += *count;
            }
        });
        KnownListConsensus {
            observed_identifier_to_best_matches: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total_count,

        }
    }

    pub fn create_balanced_bins(&self, number_of_splits: u64) -> KnownListBinSplit {
        let approximate_read_count = self.total_reads / (number_of_splits + 1);

        let mut current_total = 0;

        let mut best_match_to_original: HashMap<Vec<FastaBase>, Vec<FastaBase>> = HashMap::new();
        let mut original_to_best_match: HashMap<Vec<FastaBase>, Vec<FastaBase>> = HashMap::new();
        let mut list_to_container: HashMap<Vec<FastaBase>, usize> = HashMap::new();

        let mut container_number = 0;
        let mut dropped_read = 0;

        self.observed_identifiers_counts.iter().for_each(|(k, count)| {
            let best_hit = self.observed_identifier_to_best_matches.get(&k.clone()).unwrap();

            if best_hit.hits.len() > 0 {
                best_match_to_original.insert(best_hit.hits[0].clone(), k.clone());
                original_to_best_match.insert(k.clone(), best_hit.hits[0].clone());
                list_to_container.insert(k.clone(), container_number);

                current_total += count;
            } else {
                dropped_read += 1;
            }
            if current_total >= approximate_read_count {
                current_total = 0;
                container_number += 1;
            }
        });
        trace!("Dropped read count: {}", dropped_read);

        KnownListBinSplit { best_match_to_original, original_to_best_match, hit_to_container_number: list_to_container, bins: container_number + 1 }
    }
}

pub struct KnownListBinSplit {
    pub best_match_to_original: HashMap<Vec<FastaBase>, Vec<FastaBase>>,
    pub original_to_best_match: HashMap<Vec<FastaBase>, Vec<FastaBase>>,
    pub hit_to_container_number: HashMap<Vec<FastaBase>, usize>,
    pub bins: usize,
}

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

        let mut raw_reader = BufReader::new(File::open(filename).unwrap());

        let mut cnt = 0;

        // sort the barcodes as we go with a btree
        let mut btree = BTreeSet::new();

        for line in raw_reader.lines() {
            let bytes = line.unwrap();
            cnt += 1;
            if bytes.contains("-") {
                btree.insert(FastaBase::from_vec_u8(&bytes.split("-").collect::<Vec<&str>>().get(0).unwrap().as_bytes().to_vec()));
            } else {
                btree.insert(FastaBase::from_string(&bytes));
            }
        }
        btree
    }
}