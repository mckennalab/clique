use std::collections::{BTreeSet, HashMap, VecDeque};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::slice::Iter;
use std::sync::{Arc, Mutex};

use log::{info, warn};
use rayon::iter::ParallelBridge;
use rayon::prelude::*;
use rayon::prelude::IntoParallelRefMutIterator;

use crate::fasta_comparisons::*;
use crate::read_strategies::sequence_file_containers::*;
use crate::read_strategies::sequence_layout::*;
use crate::sorters::*;
use crate::sorters::sort_streams::SortStream;
use crate::sorters::sorter::SortStructure;
use crate::umis::sequence_clustering::*;
use crate::umis::sequence_clustering::BestHits;
use crate::utils::base_utils::edit_distance;
use crate::utils::file_utils::get_reader;
use crate::RunSpecifications;
use indicatif::ProgressBar;
use crate::sorters::balanced_split_output::RoundRobinDiskWriter;
use std::str::FromStr;
use crate::sorters::sorter::SortStructure::KNOWN_LIST;
use std::num::ParseIntError;


pub struct KnownListDiskStream {
    sorted_bins: SuperClusterOnDiskIterator,
    pattern: ReadPattern,
}

impl KnownListDiskStream {
    pub fn sort_disk_in_place(umi_type: UMIType, bin: &ReadFileContainer, sort_structure: &SortStructure, layout: &LayoutType) {
        let read_iterator = ReadIterator::new_from_bundle(bin);

        let mut sorted = Vec::new();

        for rd in read_iterator {
            let transformed_reads = transform(rd, layout);
            assert!(transformed_reads.has_original_reads());

            let field = transformed_reads.umis().unwrap().get(&umi_type).unwrap().clone();

            sorted.push((field, transformed_reads.original_reads().unwrap()));
        }
        trace!("sorted size {}", sorted.len());
        sorted.sort_by(|a, b| b.0.cmp(&a.0));

        let mut output_container = OutputReadSetWriter::from_read_file_container(&bin);
        for (string_name, read) in sorted {
            output_container.write(&read);
        }
    }
}


impl<'z> SortStream<'z> for KnownListDiskStream {
    fn from_read_iterator(read_iter: Box<dyn Iterator<Item=ReadSetContainer>>,
                          read_pattern: &ReadPattern,
                          sort_structure: &SortStructure,
                          layout: &LayoutType,
                          run_specs: &'z mut RunSpecifications) -> Self {

        match sort_structure {

            SortStructure::KNOWN_LIST { umi_type, layout_type, max_distance, on_disk, known_list } => {
                let mut consensus_manager = KnownListConsensus::new();

                let mut rrs = RoundRobinDiskWriter::from(run_specs, &read_pattern, 100);

                let bar2 = ProgressBar::new((run_specs.estimated_reads as u64));

                for rd in read_iter {
                    //trace!("Read in a read patterned {} with reads true,{},{},{}",read_pattern, rd.read_two.is_some(),rd.index_one.is_some(),rd.index_two.is_some());

                    let mut transformed_reads = transform(rd.clone(), &layout);
                    let sequence = transformed_reads.umis().unwrap().get(umi_type).unwrap().clone();

                    let corrected_hits = correct_to_known_list(&sequence, &mut known_list.lock().as_mut().unwrap(), 1);
                    consensus_manager.add_hit(&sequence, corrected_hits.clone());
                    if corrected_hits.hits.len() == 1 {
                        transformed_reads.correct_known_sequence(umi_type.clone(), corrected_hits.hits.get(0).unwrap());
                        rrs.write(&sequence, &rd);
                    }

                    bar2.inc(1);
                }
                let bins = rrs.close_and_recover_iterator();

                KnownListDiskStream { sorted_bins:  bins, pattern: read_pattern.clone() }
            }
            _ => {
                panic!("Called KnownListDiskStream using a sort structure that isn't KNOWN_LIST");
            }
        }
    }

    fn from_read_collection(read_collection: ClusteredReads, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern, run_specs: &'z mut RunSpecifications) -> Self {
        KnownListDiskStream::from_read_iterator( read_collection.reads, &pattern, sort_structure, layout, run_specs)
    }

    fn sorted_read_set(self) -> Option<SuperClusterOnDiskIterator> {
        Some(self.sorted_bins)
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
        trace!("Dropped read count: {}", dropped_read);

        KnownListBinSplit { best_match_to_original, original_to_best_match, hit_to_container_number: list_to_container, bins: container_number + 1 }
    }
}

pub struct KnownListBinSplit {
    pub best_match_to_original: HashMap<Vec<u8>, Vec<u8>>,
    pub original_to_best_match: HashMap<Vec<u8>, Vec<u8>>,
    pub hit_to_container_number: HashMap<Vec<u8>, usize>,
    pub bins: usize,
}

pub struct KnownList {
    pub name: UMIType,
    pub known_list_map: HashMap<Vec<u8>, BestHits>,
    pub known_list_subset: HashMap<Vec<u8>, Vec<Vec<u8>>>,
    pub known_list_subset_key_size: usize,
}

fn validate_barcode(barcode: &Vec<u8>) -> bool {
    barcode.iter().filter(|b| !KNOWNBASES.contains_key(*b)).map(|n| n).collect::<Vec<_>>().len() == 0 as usize
}

impl KnownList {
    const KNOWN_LIST_SEPARATOR: &'static str = ",";
    const KNOWN_TAG_SEPARATOR: &'static str = ":";

    pub fn new(knownlist_tag_and_file: &String) -> HashMap<UMIType,UMIInstance> {
        let knownlist_tag_and_files = knownlist_tag_and_file.split(KnownList::KNOWN_LIST_SEPARATOR).collect::<Vec<&str>>();

        let mut known_lists : HashMap<UMIType,UMIInstance> = HashMap::new();

        for known_list_pair in knownlist_tag_and_files {

            let knownlist_tag_and_file_vec = known_list_pair.split(KnownList::KNOWN_TAG_SEPARATOR).collect::<Vec<&str>>();

            assert_eq!(knownlist_tag_and_file_vec.len(), 3, "Unable to separate {} from {} into two tokens using divider {}", &known_list_pair, knownlist_tag_and_file, KnownList::KNOWN_TAG_SEPARATOR);

            let tag = match UMIType::from_str(knownlist_tag_and_file_vec[0]) {
                Ok(x) => {x}
                Err(_) => {panic!("Unable to parse your known list with tag: {}",&known_list_pair)}
            };

            let length = match usize::from_str(knownlist_tag_and_file_vec[1]) {
                Ok(x) => {x}
                Err(_) => {panic!("Unable to parse out UMI length from string: {}", knownlist_tag_and_file_vec[1])}
            };

            let known_list = KnownList::read_known_list_file(&tag, &knownlist_tag_and_file_vec[2], &length);

            known_lists.insert(tag.clone(),
                               UMIInstance{
                                   umi_type: tag.clone(),
                                   umi_length: length,
                                   umi_file: Some(knownlist_tag_and_file_vec[2].to_string()),
                                   known_list: Some(Arc::new(Mutex::new(known_list))),
                               });

        }
        known_lists
    }

    pub fn read_known_list_file(umi_type: &UMIType, filename: &str, starting_nmer_size: &usize) -> KnownList {

        info!("Setting up known list reader for {}", &filename);
        let btree = KnownList::create_b_tree(filename);

        // now that the barcodes are read in and in-order, we make a quick lookup prefix lookup to speed up matching later on
        let mut prefix: Option<Vec<u8>> = None;
        let mut container: Vec<Vec<u8>> = Vec::new();

        let mut existing_mapping = HashMap::new();
        let mut known_list_subset: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();

        for bytes in &btree {
            if !prefix.is_some() { prefix = Some(bytes[0..*starting_nmer_size].to_vec()); }

            existing_mapping.insert(bytes.clone(), BestHits { hits: vec![bytes.clone()], distance: 0 });

            let first_x = bytes.clone()[0..*starting_nmer_size].to_vec();
            if edit_distance(&first_x, prefix.as_ref().unwrap()) > 0 {
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

    pub fn create_b_tree(filename: &str) -> BTreeSet<Vec<u8>> {

        let mut raw_reader = get_reader(filename).unwrap();

        let mut cnt = 0;

        // sort the barcodes as we go with a btree
        let mut btree = BTreeSet::new();

        for line in raw_reader.lines() {
            let bytes = line.unwrap().as_bytes().to_vec();
            if validate_barcode(&bytes) {
                btree.insert(bytes);
            }
        }
        btree
    }

}
