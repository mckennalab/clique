use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::read_strategies::sequence_structures::*;
use crate::read_strategies::sequence_layout::*;
use crate::umis::sequence_clustering::BestHits;
use crate::utils::custom_read_sort;
use crate::sorters::sorter::ReadSortingOnDiskContainer;
use crate::sorters::sorter::SortingOutputContainer;
use crate::sorters::sorter::SortedInputContainer;
use crate::umis::sequence_clustering::*;



pub struct KnownListConsensus {
    observed_identifier_to_best_matches: HashMap<Vec<u8>, BestHits>,
    observed_identifiers_counts: HashMap<Vec<u8>, u64>,
    total_reads: u64,
}

pub struct KnownListSort {
    pub best_match_to_original: HashMap<Vec<u8>, Vec<u8>>,
    pub original_to_best_match: HashMap<Vec<u8>, Vec<u8>>,
    pub hit_to_container_number: HashMap<Vec<u8>, usize>,
    pub bins: usize,
}

impl KnownListConsensus {
    pub fn new() -> KnownListConsensus {
        let ids: HashMap<Vec<u8>, BestHits> = HashMap::new();
        let counts: HashMap<Vec<u8>, u64> = HashMap::new();
        KnownListConsensus { observed_identifier_to_best_matches: ids, observed_identifiers_counts: counts, total_reads: 0 }
    }

    pub fn add_hit(&mut self, barcode: &Vec<u8>, best_hit: BestHits) {
        self.total_reads = self.total_reads + 1;
        let new_total = *self.observed_identifiers_counts.entry(barcode.clone()).or_insert(0) + 1 as u64;
        self.observed_identifiers_counts.insert(barcode.clone(), new_total);
        self.observed_identifier_to_best_matches.insert(barcode.clone(), best_hit);
    }


    pub fn match_to_knownlist(&self) -> KnownListConsensus {
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
        KnownListConsensus {
            observed_identifier_to_best_matches: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total,
        }
    }

    pub fn filter_to_minimum_coverage(&self, minimum_coverage: u64) -> KnownListConsensus {
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
        KnownListConsensus {
            observed_identifier_to_best_matches: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total_count,

        }
    }

    pub fn read_balanced_lists(&self, number_of_splits: u64) -> KnownListSort {
        let approximate_read_count = self.total_reads / (number_of_splits + 1);

        let mut current_total = 0;

        let mut best_match_to_original: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut original_to_best_match: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut list_to_container: HashMap<Vec<u8>, usize> = HashMap::new();

        let mut container_number = 0;
        self.observed_identifiers_counts.iter().for_each(|(k, count)| {
            best_match_to_original.insert(self.observed_identifier_to_best_matches.get(&k.clone()).unwrap().hits[0].clone(), k.clone());
            original_to_best_match.insert(k.clone(), self.observed_identifier_to_best_matches.get(&k.clone()).unwrap().hits[0].clone());
            list_to_container.insert(k.clone(), container_number);

            current_total += count;

            if current_total >= approximate_read_count {
                current_total = 0;
                container_number += 1;
            }
        });
        KnownListSort { best_match_to_original, original_to_best_match, hit_to_container_number: list_to_container, bins: container_number }
    }

    pub fn consensus(read_layout: &ReadSetContainer,
                     read_files: ReadFileContainer,
                     output_file: &Path,
                     splits: &KnownListSort,
                     layout: &LayoutType) {

        assert!(layout.has_tags()); // otherwise we shouldn't be here

        let read_iterator = ReadIterator::new_from_bundle(&read_files);
        let temp_location_base_full = tempfile::tempdir().unwrap();
        let temp_location_base = temp_location_base_full.path();

        let mut temp_files = ReadSortingOnDiskContainer::create_x_bins(read_layout, splits.bins, &temp_location_base);
        let mut output_bins = temp_files.iter().map(|x| SortingOutputContainer::from_sorting_container(&x)).collect::<Vec<SortingOutputContainer>>();


        for rd in read_iterator {
            let transformed_reads = transform(rd, layout);
            assert!(transformed_reads.has_original_reads());

            let first_barcode = transformed_reads.cell_id().unwrap();

            let target_bin = splits.hit_to_container_number.get(first_barcode).unwrap();

            let original_reads = transformed_reads.original_reads().unwrap();

            output_bins[*target_bin].write_reads(original_reads);
        }


        // sort each output bin and write to the output file
        temp_files.iter().enumerate().map(|(index,bin)|{
            let sortingInput = SortedInputContainer::from_sorting_container(bin, index, &splits, &layout, &temp_location_base.to_path_buf());

        });

        //output
    }

    pub fn merge_by_umi(sic: &Vec<Box<dyn SequenceLayout>>, max_dist: u64) -> SortedInputContainer {

        let mut umi_to_reads: HashMap<Vec<u8>,Vec<ReadSetContainer>> = HashMap::new();
        sic.into_iter().for_each(|sl| {
            let umi_seq = sl.umi().clone().unwrap();
            let empty_vec: Vec<ReadSetContainer> = Vec::new();
            let mut umi_bundle = umi_to_reads.get(umi_seq).unwrap_or(&empty_vec).clone();
            umi_bundle.push(sl.original_reads().unwrap());
            umi_to_reads.insert(umi_seq.clone(), umi_bundle.clone().to_owned());
        });

        let input_list = InputList{ strings: umi_to_reads.keys().map(|k| k.clone()).collect::<Vec<Vec<u8>>>(), max_dist: 3 };
        let graph = input_list_to_graph(&input_list, string_distance, false);
        let cc_list = get_connected_components(&graph);

        let mut clusters: Vec<Vec<Vec<u8>>> = Vec::new();

        cc_list.into_iter().for_each(|cc| {

            let mut subgraph = input_list_to_graph(&InputList { strings: cc.clone(), max_dist }, string_distance, false);

            let over_set_dist = max_set_distance(&subgraph.string_to_node.keys().map(|k| k.clone()).collect::<Vec<Vec<u8>>>()) <= max_dist * 2;

            if !over_set_dist {
                clusters.push(cc);

            } else {
                let new_groups = split_subgroup(&mut subgraph);
                if new_groups.is_some() {
                    clusters.extend(new_groups.unwrap());
                } else {
                    clusters.push(cc);
                }
            }
        });
        let sorted_records: Vec<(String,Box<dyn SequenceLayout>)> = Vec::new();
        umi_to_reads.iter().for_each(|(umi, reads)| {

        });

        SortedInputContainer{ sorted_records }
    }

    pub fn consensus_from_sorted_container_10x(sic: &SortedInputContainer) -> Vec<ReadSetContainer> {
        let mut return_reads = Vec::new();
        let mut reads_to_merge = Vec::new();
        let mut current_key : Option<String> = None;

        sic.sorted_records.iter().for_each(|(barcode_string, reads)| {
            // fix this
            if !current_key.is_some() || (current_key.is_some() && string_distance(&current_key.as_ref().unwrap().clone().into_bytes(), &barcode_string.clone().into_bytes()) != 0) {
                // make a consensus for all defined reads
                if reads_to_merge.len() > 0 {

                }

            } else {
                reads_to_merge.push(reads);
            }
        });

        return_reads
    }
}
