use std::collections::{HashMap, BTreeSet, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::read_strategies::sequence_file_containers::*;
use crate::read_strategies::sequence_layout::*;
use crate::umis::sequence_clustering::BestHits;
use crate::sorters::*;
use crate::umis::sequence_clustering::*;
use crate::consensus::consensus_builders::create_seq_layout_poa_consensus;
use rayon::iter::ParallelBridge;
use rayon::prelude::IntoParallelRefMutIterator;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use log::{info, warn};
use crate::utils::file_utils::get_reader;
use crate::utils::base_utils::edit_distance;
use crate::fasta_comparisons::*;
use std::io::BufRead;
use crate::sorters::sorter::SortStructure;
use crate::sorters::sort_streams::SortStream;
use std::slice::Iter;

pub struct KnownListDiskStream {
    sorted_bins: VecDeque<PathBuf>,
    pattern: ReadPattern,
}

impl KnownListDiskStream {
    pub fn sort_disk_in_place(bin: &ReadFileContainer, sort_structure: &SortStructure, layout: &LayoutType) {

        let read_iterator = ReadIterator::new_from_bundle(bin);

        let mut sorted = Vec::new();

        for rd in read_iterator {
            let transformed_reads = transform(rd, layout);
            assert!(transformed_reads.has_original_reads());

            let field = sort_structure.get_field(&transformed_reads).unwrap();

            sorted.push((field, transformed_reads.original_reads().unwrap()));
        }
        println!("sorted size {}",sorted.len());
        sorted.sort_by(|a, b| b.0.cmp(&a.0));

        //drop(read_iterator);
        let mut output_container = OutputReadSetWriter::from_read_file_container(&bin);
        for (string_name,read) in sorted {
            output_container.write(&read);
        }
    }
}


impl SortStream for KnownListDiskStream {

    fn from_read_iterator(read_iter: ReadIterator, sort_structure: &SortStructure, layout: &LayoutType) -> Self {
        match sort_structure {
            SortStructure::KNOWN_LIST { layout_type, maximum_distance, on_disk, known_list } => {

                let mut consensus_manager = KnownListConsensus::new();
                let pattern = ReadPattern::from_read_iterator(&read_iter);
                let read_iter2 = read_iter.new_reset();
                let read_iter3 = read_iter.new_reset();
                println!("round 1");
                for rd in read_iter {
                    let transformed_reads = transform(rd, &layout);
                    let sequence = sort_structure.get_field(&transformed_reads).unwrap();

                    let corrected_hits = correct_to_known_list(&sequence, &mut known_list.lock().as_mut().unwrap(), 1);
                    consensus_manager.add_hit(&sequence, corrected_hits);
                }

                let kcl = consensus_manager.match_to_knownlist();
                let splits = consensus_manager.create_balanced_bins(10 as u64);

                let temp_location_base_full = tempfile::tempdir().unwrap();
                let temp_location_base = Path::new("./tmp/");

                let mut temp_files = OutputReadSetWriter::create_x_bins(&read_iter3, &"unsorted".to_string(), splits.bins, &temp_location_base);
                let mut output_bins = temp_files.iter().map(|(id,x)| {
                    println!("opening {}",x.read_one.to_str().unwrap());
                    OutputReadSetWriter::from_read_file_container(&x)}).collect::<Vec<OutputReadSetWriter>>();

                let mut counts: HashMap<usize,i32> = HashMap::new();
                println!("round 2");
                let mut sorted_reads = 0;
                for rd in read_iter2 {
                    let transformed_reads = transform(rd, &layout);
                    let sequence = sort_structure.get_field(&transformed_reads).unwrap();

                    let target_bin = splits.hit_to_container_number.get(&sequence);

                    if let Some(target_bin) = target_bin {
                        sorted_reads += 1;
                        let original_reads = transformed_reads.original_reads().unwrap();
                        let bin_count: i32 = *counts.get(target_bin).unwrap_or(&0);
                        counts.insert(*target_bin,bin_count + 1);
                        output_bins[*target_bin].write(&original_reads);
                    }
                }
                output_bins.iter().for_each(|x| x.print_read_count());
                for x in counts.iter() {
                    println!("Bin {} wrote {}",x.0, x.1);
                }
                // make sure we flush all the buffers
                drop(output_bins);
                println!("sort files {}",sorted_reads);
                temp_files.iter().for_each(|(size,temp_file)| {
                    KnownListDiskStream::sort_disk_in_place(&temp_file, sort_structure, layout);
                });

                let bins = VecDeque::from(temp_files.iter().map(|n| n.1.read_one.clone()).collect::<Vec<PathBuf>>());
                println!("Bin size: {}", &bins.len());
                KnownListDiskStream{ sorted_bins: bins, pattern }
            }
            _ => {
                panic!("Called KnownListDiskStream using a sort structure that isn't KNOWN_LIST");
            }
        }


    }

    fn from_read_collection(read_collection: ReadCollection, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern) -> Self {
        KnownListDiskStream::from_read_iterator(ReadIterator::from_collection(read_collection),sort_structure,layout)
    }

    fn sorted_read_set(&mut self) -> Option<ReadCollectionIterator> {
        match self.sorted_bins.len() {
            0 => None,
            _ => {
                let next_file = self.sorted_bins.pop_front();
                match next_file {
                    None => None,
                    Some(x) => Some(ReadCollectionIterator::new_from_file(vec![x], self.pattern.clone())),
                }
            }
        }
    }
}

pub struct KnownListConsensus {
    observed_identifier_to_best_matches: HashMap<Vec<u8>, BestHits>,
    observed_identifiers_counts: HashMap<Vec<u8>, u64>,
    total_reads: u64,
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

        info!("total = {}, unresolved = {}, resolved = {}, final = {}", total, unresolved, resolved, clean_mapping.len());

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

    pub fn create_balanced_bins(&self, number_of_splits: u64) -> KnownListBinSplit {
        let approximate_read_count = self.total_reads / (number_of_splits + 1);

        let mut current_total = 0;

        let mut best_match_to_original: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut original_to_best_match: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut list_to_container: HashMap<Vec<u8>, usize> = HashMap::new();

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
        println!("Dropped read count: {}",dropped_read);

        KnownListBinSplit { best_match_to_original, original_to_best_match, hit_to_container_number: list_to_container, bins: container_number + 1}
    }
}

pub struct KnownListBinSplit {
    pub best_match_to_original: HashMap<Vec<u8>, Vec<u8>>,
    pub original_to_best_match: HashMap<Vec<u8>, Vec<u8>>,
    pub hit_to_container_number: HashMap<Vec<u8>, usize>,
    pub bins: usize,
}

pub struct KnownList {
    pub name: String,
    pub known_list_map: HashMap<Vec<u8>, BestHits>,
    pub known_list_subset: HashMap<Vec<u8>, Vec<Vec<u8>>>,
    pub known_list_subset_key_size: usize,
}

fn validate_barcode(barcode: &Vec<u8>) -> bool {
    barcode.iter().filter(|b| !KNOWNBASES.contains_key(*b)).map(|n| n).collect::<Vec<_>>().len() == 0 as usize
}

impl KnownList {
    const KNOWN_LIST_SEPARATOR: &'static str = ":";

    pub fn new(knownlist_tag_and_file: &String, starting_nmer_size: usize) -> KnownList {

        // TODO: fix the sep here
        let knownlist_tag_and_file_vec = knownlist_tag_and_file.split(":").collect::<Vec<&str>>();
        assert_eq!(knownlist_tag_and_file_vec.len(), 2);

        let mut existing_mapping = HashMap::new();
        let mut known_list_subset: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();

        let mut raw_reader = get_reader(knownlist_tag_and_file_vec[1]).unwrap();

        let mut cnt = 0;

        let mut btree = BTreeSet::new();
        println!("Adding known barcodes...");
        for line in raw_reader.lines() {
            let bytes = line.unwrap().as_bytes().to_vec();
            if validate_barcode(&bytes) {
                btree.insert(bytes);
            }
        }

        // now the barcodes are in order, use this to generate grouped keys
        let mut prefix: Option<Vec<u8>> = None;
        let mut container: Vec<Vec<u8>> = Vec::new();
        for bytes in &btree {
            if !prefix.is_some() { prefix = Some(bytes[0..starting_nmer_size].to_vec()); }

            existing_mapping.insert(bytes.clone(), BestHits { hits: vec![bytes.clone()], distance: 0 });

            let first_x = bytes.clone()[0..starting_nmer_size].to_vec();
            if edit_distance(&first_x, prefix.as_ref().unwrap()) > 0 {
                known_list_subset.insert(prefix.unwrap(), container.clone());
                container.clear();
                prefix = Some(first_x);
            }
            container.push(bytes.clone());
        }
        known_list_subset.insert(prefix.unwrap(), container.clone());
        KnownList {
            name: knownlist_tag_and_file_vec[0].to_string(),
            known_list_map: existing_mapping,
            known_list_subset,
            known_list_subset_key_size: starting_nmer_size,
        }
    }
}
