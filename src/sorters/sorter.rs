use std::borrow::BorrowMut;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::slice::Iter;
use std::sync::{Arc, Mutex};

use flate2::Compression;
use flate2::write::GzEncoder;
use rust_htslib::bgzf::Writer;

use log::{info, warn};

use crate::consensus::consensus_builders::create_seq_layout_poa_consensus;
use crate::read_strategies::sequence_file_containers::{ReadFileContainer, ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::OutputReadSetWriter;
use crate::read_strategies::sequence_file_containers::ReadCollection;
use crate::read_strategies::sequence_layout::LayoutType;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::read_strategies::sequence_layout::transform;
use crate::sorters::known_list::KnownList;
use crate::sorters::known_list::KnownListBinSplit;
use crate::sorters::known_list::KnownListConsensus;
use crate::sorters::known_list::KnownListDiskStream;
use crate::sorters::sort_streams::*;
use crate::read_strategies::sequence_file_containers::ReadCollectionIterator;


#[derive(Clone)]
pub enum SortStructure {
    KNOWN_LIST { layout_type: LayoutType, maximum_distance: usize, on_disk: bool, known_list: Arc<Mutex<KnownList>> },
    HD_UMI { layout_type: LayoutType, on_disk: bool },
    LD_UMI { layout_type: LayoutType, on_disk: bool },
}

impl SortStructure {
    pub fn from_layout(layout: &LayoutType, known_lists: HashMap<LayoutType, Arc<Mutex<KnownList>>>) -> Vec<SortStructure> {
        match layout {
            LayoutType::TENXV3 => {
                let mut ret = Vec::new();
                ret.push(SortStructure::KNOWN_LIST { layout_type: LayoutType::TENXV3, maximum_distance: 1, on_disk: true, known_list: known_lists.get(&LayoutType::TENXV3).unwrap().clone() });
                ret.push(SortStructure::LD_UMI { layout_type: LayoutType::TENXV3, on_disk: false });
                ret
            }
            LayoutType::TENXV2 => {
                let mut ret = Vec::new();
                ret.push(SortStructure::KNOWN_LIST { layout_type: LayoutType::TENXV2, maximum_distance: 1, on_disk: true, known_list: known_lists.get(&LayoutType::TENXV2).unwrap().clone() });
                ret.push(SortStructure::LD_UMI { layout_type: LayoutType::TENXV2, on_disk: false });
                ret
            }
            LayoutType::PAIREDUMI => { unimplemented!() }
            LayoutType::PAIRED => { unimplemented!() }
            LayoutType::SCI => { unimplemented!() }
            _ => { unimplemented!() }
        }
    }

    /// I couldn't get my sort structure to take a generic 'function parameter' to call for the correct sequences, so this
    /// method has to do that work in a less elegant way. Maybe I'll figure it out later
    ///
    pub fn get_known_list_value(&self, layout_type: &LayoutType, seq_layout: &dyn SequenceLayout, known_lists: HashMap<LayoutType, PathBuf>) -> PathBuf {
        match self {
            SortStructure::KNOWN_LIST { layout_type, maximum_distance, on_disk, known_list } => {
                match layout_type {
                    LayoutType::TENXV3 => { known_lists.get(layout_type).unwrap().clone() }
                    LayoutType::TENXV2 => { unimplemented!() }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
            SortStructure::HD_UMI { layout_type, on_disk } => {
                unimplemented!()
            }
            SortStructure::HD_UMI { layout_type, on_disk } => {
                unimplemented!()
            }
            _ => { unimplemented!() }
        }
    }

    pub fn get_field(&self, seq_layout: &impl SequenceLayout) -> Option<Vec<u8>> {
        match self {
            SortStructure::KNOWN_LIST { layout_type, maximum_distance, on_disk, known_list } => {
                match layout_type {
                    LayoutType::TENXV3 => { Some(seq_layout.cell_id().unwrap().clone()) }
                    LayoutType::TENXV2 => { Some(seq_layout.cell_id().unwrap().clone()) }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
            SortStructure::HD_UMI { layout_type, on_disk } => {
                match layout_type {
                    LayoutType::TENXV3 => { Some(seq_layout.umi().unwrap().clone()) }
                    LayoutType::TENXV2 => { Some(seq_layout.umi().unwrap().clone()) }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
            SortStructure::LD_UMI { layout_type, on_disk } => {
                match layout_type {
                    LayoutType::TENXV3 => { Some(seq_layout.umi().unwrap().clone()) }
                    LayoutType::TENXV2 => { Some(seq_layout.umi().unwrap().clone()) }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
        }
    }
}

pub struct Sorter {}

impl Sorter {
    pub fn sort(sort_list: Vec<SortStructure>, input_reads: &ReadFileContainer, tmp_location: &String, sorted_output: &String, layout: &LayoutType) -> Vec<ReadIterator> {
        let temp_location_base = Path::new(tmp_location);
        println!("Sorting reads1...");

        let mut read_iterator = ReadIterator::new_from_bundle(input_reads);

        let mut current_iterators = Vec::new();

        current_iterators.push(read_iterator);

        println!("Sorting reads...");
        for sort in sort_list {
            let mut next_level_iterators = Vec::new();

            for mut iter in current_iterators {
                match Sorter::sort_level(&sort, iter, layout) {
                    None => {}
                    Some(x) => {next_level_iterators.extend(x);}
                }
            }

            current_iterators = next_level_iterators;
            println!("Done sorting reads...length {}",current_iterators.len());

        }
        println!("Done sorting reads...length {}",current_iterators.len());
        current_iterators
    }

    pub fn sort_level(sort_structure: &SortStructure, iterator: ReadIterator, layout: &LayoutType) -> Option<Vec<ReadIterator>> {
        match sort_structure {
            SortStructure::KNOWN_LIST { layout_type, maximum_distance, on_disk, known_list } => {
                assert_eq!(*on_disk,true);
                let mut sorter = KnownListDiskStream::from_read_iterator(iterator, sort_structure, layout);
                let mut read_sets = sorter.sorted_read_set();
                match read_sets {
                    None => None,
                    Some(x) => {
                        let mut ret = Vec::new();
                        for ci in x {
                            ret.push(ReadIterator::from_collection(ci));
                        }
                        Some(ret)
                    }
                }
            }
            SortStructure::HD_UMI { layout_type, on_disk} |
            SortStructure::LD_UMI { layout_type, on_disk } => {
                assert_eq!(*on_disk,false);
                let mut sorter = ClusteredMemorySortStream::from_read_iterator(iterator, sort_structure, layout);
                let read_sets = sorter.sorted_read_set();
                match read_sets {
                    None => None,
                    Some(x) => {
                        let mut ret = Vec::new();
                        for ci in x {
                            ret.push(ReadIterator::from_collection(ci));
                        }
                        Some(ret)
                    }
                }
            }
        }


    }
    /*
        pub fn bin_by_tag(bin: &ReadFileContainer, sort_structure: &SortStructure, layout: &LayoutType) -> Vec<(usize, ReadFileContainer)> {
            let read_iterator = ReadIterator::new_from_bundle(bin);

            let temp_location_base = Path::new("./tmp/");
            let mut temp_files = OutputReadSetWriter::create_x_bins(&read_iterator, &"unsorted".to_string(), splits.bins, &temp_location_base);
            let mut output_bins = temp_files.iter().map(|(id,x)|
                OutputReadSetWriter::from_sorting_container(&x)).collect::<Vec<OutputReadSetWriter>>();

            let mut counts: HashMap<usize,i32> = HashMap::new();

            let mut cnt = 0;
            let mut dropped = 0;
            for rd in read_iterator {
                let transformed_reads = transform(rd, layout);
                assert!(transformed_reads.has_original_reads());

                let field = sort_structure.get_field(&transformed_reads).unwrap();

                let target_bin = splits.hit_to_container_number.get(field);

                if let Some(target_bin) = target_bin {
                    let original_reads = transformed_reads.original_reads().unwrap();
                    cnt += 1;
                    let bin_count: i32 = *counts.get(target_bin).unwrap_or(&0);
                    counts.insert(*target_bin,bin_count + 1);
                    output_bins[*target_bin].write(original_reads);
                } else {
                    dropped += 1;
                }
            }
            counts.iter().for_each(|(k,v)| println!("K {} v {}",k,v));

            info!("Read count {}, dropped {}",cnt, dropped);

            temp_files
        }
    */
}
