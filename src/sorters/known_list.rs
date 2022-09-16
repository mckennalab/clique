use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::read_strategies::sequence_structures::*;
use crate::read_strategies::sequence_layout::*;
use crate::umis::sequence_clustering::BestHits;
use crate::utils::custom_read_sort;
use crate::sorters::sorter::ReadSortingOnDiskContainer;
use crate::sorters::sorter::SortingOutputContainer;
use crate::sorters::sorter::SortedInputContainer;
use crate::umis::sequence_clustering::*;
use crate::consensus::consensus_builders::create_seq_layout_poa_consensus;
use rayon::iter::ParallelBridge;
use rayon::prelude::IntoParallelRefMutIterator;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

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
            let best_hit = self.observed_identifier_to_best_matches.get(&k.clone()).unwrap();
            if best_hit.hits.len() > 0 {
                best_match_to_original.insert(best_hit.hits[0].clone(), k.clone());
                original_to_best_match.insert(k.clone(), best_hit.hits[0].clone());
                list_to_container.insert(k.clone(), container_number);

                current_total += count;
            } else {
                println!("Missing best hit for {} hits {:?}",String::from_utf8(k.clone()).unwrap(),best_hit.hits);
            }
            if current_total >= approximate_read_count {
                current_total = 0;
                container_number += 1;
            }
        });
        KnownListSort { best_match_to_original, original_to_best_match, hit_to_container_number: list_to_container, bins: container_number + 1}
    }

    pub fn consensus(read_file_bundle: &ReadFileContainer,
                     splits: &KnownListSort,
                     layout: &LayoutType,
                     output_file: &PathBuf,
                     threads: usize) {

        assert!(layout.has_tags()); // otherwise we shouldn't be here

        let temp_location_base_full = tempfile::tempdir().unwrap();
        let temp_location_base = Path::new("./tmp/");
        let read_iterator = ReadIterator::new_from_bundle(read_file_bundle);
        let mut temp_files = ReadSortingOnDiskContainer::create_x_bins(&read_iterator, splits.bins, &temp_location_base);
        let mut output_bins = temp_files.iter().map(|(id,x)| SortingOutputContainer::from_sorting_container(&x)).collect::<Vec<SortingOutputContainer>>();

        let mut counts: HashMap<usize,i32> = HashMap::new();

        let mut cnt = 0;
        let mut dropped = 0;
        for rd in read_iterator {
            let transformed_reads = transform(rd, layout);
            assert!(transformed_reads.has_original_reads());

            let first_barcode = transformed_reads.cell_id().unwrap();

            let target_bin = splits.hit_to_container_number.get(first_barcode);

            if let Some(target_bin) = target_bin {
                let original_reads = transformed_reads.original_reads().unwrap();
                cnt += 1;
                let bin_count: i32 = *counts.get(target_bin).unwrap_or(&0);
                counts.insert(*target_bin,bin_count + 1);
                output_bins[*target_bin].write_reads(original_reads);
            } else {
                dropped += 1;
            }
        }
        counts.iter().for_each(|(k,v)| println!("K {} v {}",k,v));

        println!("COunts {} dropped {}",cnt, dropped);

        output_bins.iter_mut().for_each(|o| o.close());
        let mut output_file = BufWriter::new(File::create(output_file).unwrap());
        let output = Arc::new(Mutex::new(output_file));

        // sort each output bin and write to the output file
        //rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

        temp_files.par_iter().for_each(|(index,bin)|{
            println!("sorting input ");
            let sortingInput = SortedInputContainer::from_sorting_container(bin, *index, &splits, &layout, &temp_location_base.to_path_buf());
            println!("sorting input size {}",sortingInput.sorted_records.len());
            let consensuses = KnownListConsensus::consensus_from_sorted_container_10x(&sortingInput, 1);

            let output = Arc::clone(&output);
            let mut output = output.lock().unwrap();

            for con in consensuses {
                write!(output,"R1={}\n,R2={}\n,I1={}\n,I2={}\n",
                       String::from_utf8(con.read_one).unwrap(),
                    if con.read_two.is_some() {String::from_utf8(con.read_two.unwrap()).unwrap()} else {"NONE".to_string()},
                       if con.index_one.is_some() {String::from_utf8(con.index_one.unwrap()).unwrap()} else {"NONE".to_string()},
                       if con.index_two.is_some() {String::from_utf8(con.index_two.unwrap()).unwrap()} else {"NONE".to_string()});
            }
        });
        let output = Arc::clone(&output);
        let mut output = output.lock().unwrap();
        output.flush();
    }

    pub fn submerge_by_umi(sic: &Vec<&Box<dyn SequenceLayout>>, max_dist: u64) -> Vec<SequenceSetContainer> {

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
        let mut sorted_records= Vec::new();
        umi_to_reads.iter().for_each(|(umi, reads)| {
            let consensus = create_seq_layout_poa_consensus(reads);
            sorted_records.push(consensus);
        });
        sorted_records
    }

    pub fn consensus_from_sorted_container_10x(sic: &SortedInputContainer, max_distance: u64) -> Vec<SequenceSetContainer> {
        let mut return_reads = Vec::new();
        let mut reads_to_merge = Vec::new();
        let mut current_key : Option<String> = None;

        sic.sorted_records.iter().for_each(|(barcode_string, reads)| {
            // fix this
            if !current_key.is_some() || (current_key.is_some() && string_distance(&current_key.as_ref().unwrap().clone().into_bytes(), &barcode_string.clone().into_bytes()) != 0) {
                // make a consensus for all defined reads
                if reads_to_merge.len() > 0 {
                    let consensus = KnownListConsensus::submerge_by_umi(&reads_to_merge, max_distance);
                    return_reads.extend(consensus);
                }
                reads_to_merge.clear();
                current_key = Some(barcode_string.clone());
            }
            reads_to_merge.push(reads);

        });

        if reads_to_merge.len() > 0 {
            let consensus = KnownListConsensus::submerge_by_umi(&reads_to_merge, max_distance);
            return_reads.extend(consensus);
        }
        return_reads
    }
}
