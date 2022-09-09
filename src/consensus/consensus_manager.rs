use std::collections::HashMap;
use std::path::Path;

use utils::custom_read_sort;

use crate::read_strategies::sequence_layout::{LayoutType, ReadFileContainer, ReadIterator};
use crate::umis::sequenceclustering::BestHits;
use read_strategies::sequence_layout::transform;

pub struct ConsensusManager {
    observed_identifier_to_best_matches: HashMap<Vec<u8>, BestHits>,
    observed_identifiers_counts: HashMap<Vec<u8>, u64>,
    total_reads: u64,
}


pub struct ConsensusManagerSort {
    one_to_one_best_match_to_original_barcode: HashMap<Vec<u8>, Vec<u8>>,
    hit_to_container_number: HashMap<Vec<u8>, usize>,
}


impl ConsensusManager {
    pub fn new() -> ConsensusManager {
        let ids: HashMap<Vec<u8>, BestHits> = HashMap::new();
        let counts: HashMap<Vec<u8>, u64> = HashMap::new();
        ConsensusManager { observed_identifier_to_best_matches: ids, observed_identifiers_counts: counts, total_reads: 0 }
    }

    pub fn add_hit(&mut self, barcode: &Vec<u8>, best_hit: BestHits) {
        self.total_reads = self.total_reads + 1;
        let new_total = *self.observed_identifiers_counts.entry(barcode.clone()).or_insert(0) + 1 as u64;
        self.observed_identifiers_counts.insert(barcode.clone(), new_total);
        self.observed_identifier_to_best_matches.insert(barcode.clone(), best_hit);
    }


    pub fn unified_consensus_list(&self) -> ConsensusManager {
        let mut total = 0;
        let mut unresolved = 0;
        let mut resolved = 0;
        let mut clean_mapping: HashMap<Vec<u8>, BestHits> = HashMap::new();
        let mut clean_ids: HashMap<Vec<u8>, u64> = HashMap::new();

        self.observed_identifier_to_best_matches.iter().for_each(|(k, v)| {
            let count = *self.observed_identifiers_counts.get(&k.clone()).unwrap();
            total += count;
            if v.hits.len() > 1 {
                let matching: Vec<&Vec<u8>> = v.hits.iter().filter(|hit| self.observed_identifiers_counts.contains_key(*hit)).collect::<Vec<&Vec<u8>>>();
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
        println!("total = {}, unresolved = {}, resolved = {}, final = {}", total, unresolved, resolved, clean_mapping.len());
        ConsensusManager {
            observed_identifier_to_best_matches: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total,
        }
    }

    pub fn filter_to_minimum_coverage(&self, minimum_coverage: u64) -> ConsensusManager {
        let mut clean_mapping: HashMap<Vec<u8>, BestHits> = HashMap::new();
        let mut clean_ids: HashMap<Vec<u8>, u64> = HashMap::new();

        let mut total_count = 0;

        self.observed_identifiers_counts.iter().map(|(k, count)| {
            if *count >= minimum_coverage {
                clean_ids.insert(k.clone(), *count);
                clean_mapping.insert(k.clone(), self.observed_identifier_to_best_matches.get(k).unwrap().clone());
                total_count += *count;
            }
        });
        ConsensusManager {
            observed_identifier_to_best_matches: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total_count,
        }
    }

    pub fn create_split_manager(reverse_list: &HashMap<Vec<u8>, Vec<u8>>, list_to_container: &HashMap<Vec<u8>, usize>) -> ConsensusManagerSort {
        ConsensusManagerSort {
            one_to_one_best_match_to_original_barcode: reverse_list.to_owned(),
            hit_to_container_number: list_to_container.to_owned(),
        }
    }

    pub fn read_balanced_lists(&self, number_of_splits: u64) -> Vec<ConsensusManagerSort> {
        let approximate_read_count = self.total_reads / (number_of_splits + 1);
        let mut splits = Vec::new();

        let mut current_total = 0;

        let mut reverse_list: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut list_to_container: HashMap<Vec<u8>, usize> = HashMap::new();

        let mut container_number = 0;
        self.observed_identifiers_counts.iter().for_each(|(k, count)| {
            if !clean_mapping.is_some() {
                list_to_container = Some(HashMap::new());
                reverse_list = Some(HashMap::new());
                current_total = 0;
            }

            if let Some(x) = &mut clean_mapping {
                x.insert(k.clone(), self.observed_identifier_to_best_matches.get(&k.clone()).unwrap().clone());
            }

            if let Some(x) = &mut clean_ids {
                x.insert(k.clone(), *count);
            }

            if let Some(x) = &mut reverse_list {
                assert_eq!(self.observed_identifier_to_best_matches.get(&k.clone()).unwrap().hits.len(), 1);
                x.insert(self.observed_identifier_to_best_matches.get(&k.clone()).unwrap().hits[0].clone(), k.clone());
            }

            current_total += count;

            if current_total >= approximate_read_count {
                splits.push(ConsensusManager::create_split_manager(
                    reverse_list.as_ref().unwrap(),
                    list_to_container.as_ref().unwrap(),
                ));
                list_to_container = None;
                reverse_list = None;
                container_number += 1;
            }
        });

        if reverse_list.is_some() {
            splits.push(ConsensusManager::create_split_manager(
                reverse_list.as_ref().unwrap(),
                list_to_container.as_ref().unwrap(),
            ));
            list_to_container = None;
            reverse_list = None;
            container_number += 1;
        }

        println!("cnt = {} and {}", cnt, self.observed_identifiers_counts.len());
        splits
    }

    pub fn resort_reads(read_files: ReadFileContainer,
                        output_dir: &Path,
                        output_prefix: &String,
                        temp_location_base: &Path,
                        splits: Vec<ConsensusManagerSort>,
                        layout: &LayoutType) -> bool {

        let read_iterator = ReadIterator::new_from_bundle(&read_files);
        let tmp_files = custom_read_sort::create_tmp_list(splits.len());

        assert!(layout.has_identifying_sequence()); // otherwise we shouldn't be here

        for rd in read_iterator {

            let transformed_reads = transform(rd, &read_layout);
            let first_hit = transformed_reads.get_unique_sequences().unwrap()[0].clone();
            target_bin =

        }
        println!("Count = {}, mto = {}", cnt, mto);
        consensus_manager.unified_consensus_list();
        let splits = consensus_manager.read_balanced_lists(10);
        println!("sizes = {}", splits.len());
    }
}