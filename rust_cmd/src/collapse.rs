use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::path::{Path, PathBuf};
use std::string;
use crate::alignment_functions::{align_reads, fast_align_reads};
use crate::InstanceLivedTempDir;
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration, UMISortType};
use crate::reference::fasta_reference::ReferenceManager;
use bio::io::fasta::*;
use indicatif::ProgressBar;
use itertools::Itertools;
use log::Level::Warn;
use shardio::{DefaultSort, Range, ShardReader, ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};
use crate::sequence_lookup::KnownLookup;
use rustc_hash::FxHashMap;
use typetag::__private::inventory::iter;
use crate::consensus::consensus_builders::write_consensus_reads;
use crate::umis::known_list::KnownList;
use crate::umis::sequence_clustering::{correct_to_known_list, get_connected_components, input_list_to_graph, InputList, split_subgroup, string_distance};

pub fn collapse(reference: &String,
                final_output: &String,
                temp_directory: &mut InstanceLivedTempDir,
                read_structure: &SequenceLayoutDesign,
                max_reference_multiplier: &f64,
                min_read_length: &usize,
                read1: &String,
                read2: &String,
                index1: &String,
                index2: &String,
                threads: &usize) {


    // load up the reference files
    let rm = ReferenceManager::from(&reference, 8);

    // validate that each reference has the specified capture groups
    let validated_references = rm.references.iter().
        map(|rf| read_structure.validate_reference_sequence(&rf.1.sequence_u8)).all(|x| x == true);

    assert!(validated_references, "The reference sequences do not match the capture groups specified in the read structure file.");

    // align whatever combination of reads to the reference, collapsing down to a single sequence
    let aligned_temp = temp_directory.temp_file("aligned_reads.fasta");

    info!("Aligning reads to reference in directory {}", aligned_temp.as_path().display());

    let ret = fast_align_reads(&true,
                               read_structure,
                               &false,
                               &false,
                               &rm,
                               aligned_temp.as_path(),
                               max_reference_multiplier,
                               min_read_length,
                               read1,
                               read2,
                               index1,
                               index2,
                               threads);

    info!("Sorting the aligned reads");

    let mut read_count = ret.0;
    let mut sorted_input = ret.1;

    info!("Sorting by read tags");

    // sort the reads by the tags
    let mut levels = 0;
    read_structure.get_sorted_umi_configurations().iter().for_each(|tag| {
        match tag.sort_type {
            UMISortType::KnownTag => {
                let ret = sort_known_level(temp_directory, &sorted_input, &tag, &read_count);
                sorted_input = ret.1;
                read_count = ret.0;
            }
            UMISortType::DegenerateTag => {
                let ret = sort_degenerate_level(temp_directory, &sorted_input, &tag, &levels, &read_count);
                sorted_input = ret.1;
                read_count = ret.0;
            }
        }
        levels += 1;
    });
    info!("writing consensus reads");

    // collapse the final reads down to a single sequence and write everything to the disk
    write_consensus_reads(&sorted_input, final_output, &threads, levels, &read_count, &rm, &40);
}


fn consensus(input: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut consensus = Vec::new();

    for i in 0..input[0].len() {
        let mut counter = HashMap::new();

        input.iter().for_each(|vector| {
            *counter.entry(&vector[i]).or_insert(0) += 1;
        });

        let mut max = 0;
        let mut consensus_byte = b'N';

        for (byte, count) in counter {
            if count > max {
                max = count;
                consensus_byte = *byte;
            }
        }

        consensus.push(consensus_byte);
    }

    consensus
}

struct TagStrippingDiskBackedBuffer {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
    writen_reads: usize,
    shard_writer: Option<Box<ShardWriter<SortingReadSetContainer>>>,
    shard_sender: Option<Box<ShardSender<SortingReadSetContainer>>>,
    output_file: PathBuf,
    tag: UMIConfiguration,
    hash_map: FxHashMap<String, usize>,
}

impl TagStrippingDiskBackedBuffer {
    pub fn new(output_file: PathBuf, max_size: usize, tag: UMIConfiguration) -> TagStrippingDiskBackedBuffer {
        TagStrippingDiskBackedBuffer {
            buffer: VecDeque::new(),
            max_buffer_size: max_size,
            writen_reads: 0,
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
        }
    }

    /// Pushes a new item onto the buffer
    /// If the buffer is full, it will write the buffer to disk and clear it.
    ///
    pub fn push(&mut self, mut item: SortingReadSetContainer) {
        let key_value = item.ordered_unsorted_keys.pop_front().unwrap();
        item.ordered_unsorted_keys.push_front(key_value.clone()); // we want to keep the key in the list for now, we'll remove it later
        assert_eq!(key_value.0, self.tag.symbol);
        *self.hash_map.entry(FastaBase::to_string(&key_value.1)).or_insert(0) += 1;

        match (&self.shard_writer, self.writen_reads >= self.max_buffer_size) {
            (None, true) => {
                self.buffer.push_back(item);
                self.dump_buffer();
            }
            (None, false) => {
                self.buffer.push_back(item);
            }
            (Some(_), _) => {
                self.shard_sender.as_mut().unwrap().send(item).unwrap();
            }
        }
    }

    pub fn dump_buffer(&mut self) {
        self.shard_writer = Some(Box::new(ShardWriter::new(&self.output_file,
                                                           32,
                                                           256,
                                                           1 << 16).unwrap()));
        self.shard_sender = Some(Box::new(self.shard_writer.as_mut().unwrap().get_sender()));
        self.buffer.iter().for_each(|x| {
            self.shard_sender.as_mut().unwrap().send(x.clone()).unwrap();
        });
    }

    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        let string_set = Vec::from_iter(self.hash_map.keys()).iter().map(|s| s.as_bytes().to_vec()).collect::<Vec<Vec<u8>>>();
        let collection = InputList { strings: string_set, max_dist: u64::try_from(self.tag.max_distance).unwrap() };
        let progress_bar = self.hash_map.len() > 50000;
        let graph = input_list_to_graph(&collection, string_distance, progress_bar);

        info!("creating connected components for degenerate sequences");
        let cc = get_connected_components(&graph);
        info!("raw connected components has {} components",cc.len());
        let mut final_correction: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        for group in cc {
            let minilist = InputList { strings: group, max_dist: u64::try_from(self.tag.max_distance.clone()).unwrap() };
            /*let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let is_subgroups = split_subgroup(&mut minigraph);

            match is_subgroups {
                None => {*/
                    let group = minilist.strings.clone();
                    for s in minilist.strings {
                        final_correction.insert(s.clone(), consensus(&group));
                    };
            /*
                }
                Some(x) => {
                    for sgroup in &x {
                        for s in sgroup {
                            final_correction.insert(s.clone(), consensus(&sgroup));
                        }
                    }
                }
            }*/
        }
        final_correction
    }

    pub fn close(&mut self) -> (VecDeque<SortingReadSetContainer>, FxHashMap<Vec<u8>, Vec<u8>>) {
        let final_correction = self.correct_list();
        //std::mem::replace(&mut response.headers, hyper::header::Headers::new())
        let ret = (std::mem::replace(&mut self.buffer, VecDeque::new()), final_correction);
        self.buffer = VecDeque::new();
        self.hash_map.clear();
        ret
    }

    pub fn clean(&mut self) {
        self.buffer.clear();
        self.hash_map.clear();
    }
}


/// Sorts the reads by the known tag
///
/// we group reads into a container where previous tags all match. We then determine the clique of
/// degenerate tags within the container and correct the sequences to the consensuses of cliques within
/// the group. For example, if we've sorted by a 10X cell ID tag, that we would collect cells that have the
/// same cell ID, and then correct the UMI sequences within the cell to the consensuses of the UMI sequences
///
/// # Arguments
///     * `temp_directory` - The temporary directory to use for sorting
///    * `reader` - The reader to sort
///   * `tag` - The tag to sort by
///
/// # Returns
///    * `ShardReader` - The sorted reader
///
/// # Errors
///     * `std::io::Error` - If there is an error writing to the temporary directory
///
/// # Panics
///    * If the tag is not a known tag
///
/// # Examples
///
pub fn sort_degenerate_level(temp_directory: &mut InstanceLivedTempDir,
                             reader: &ShardReader<SortingReadSetContainer>,
                             tag: &UMIConfiguration,
                             iteration: &usize, read_count: &usize) -> (usize, ShardReader<SortingReadSetContainer>) {
    info!("Sorting degenerate level {}",tag.symbol);
    // create a new output
    let mut processed_reads = 0;
    let mut added_to_group = 0;
    let mut sent_reads = 0;
    let mut dropped_collection = 0;

    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&aligned_temp, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();

    let mut sender = sharded_output.get_sender();
    let bar = ProgressBar::new(read_count.clone() as u64);

    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut current_sorting_bin: TagStrippingDiskBackedBuffer =
        TagStrippingDiskBackedBuffer::new(temp_directory.temp_file(format!("{}_{}.fasta", 0, tag.order).as_str()), 1000000, tag.clone());

    let mut largest_group_size = 0;
    let mut group_size = 0;

    reader.iter_range(&Range::all()).unwrap().enumerate().for_each(|(i, y)| {
        bar.inc(1);
        let mut x = y.as_ref().unwrap().clone();
        let mut x2 = y.as_ref().unwrap().clone();

        processed_reads += 1;

        if !last_read.is_none() {
            assert_eq!(x.ordered_sorting_keys.len(), *iteration);
            if last_read.as_ref().unwrap().cmp(&mut x) == Ordering::Equal {
                current_sorting_bin.push(x);
                group_size += 1;
                added_to_group += 1;
            } else {
                if (tag.maximum_subsequences.is_some() && tag.maximum_subsequences.unwrap() > current_sorting_bin.buffer.len()) || tag.maximum_subsequences.is_none() {
                    let (buffer, correction) = current_sorting_bin.close();
                    if correction.len() > largest_group_size {
                        largest_group_size = correction.len();
                    }
                    for mut y in buffer {
                        let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
                        let corrected = correction.get(&FastaBase::to_vec_u8(&key_value.1)).unwrap();
                        y.ordered_sorting_keys.push((key_value.0, FastaBase::from_vec_u8(corrected)));
                        sender.send(y).unwrap();
                        sent_reads += 1;
                    }
                } else {
                    dropped_collection += 1;
                }
                group_size = 0;
                current_sorting_bin.clean();
            }
        } else {
            current_sorting_bin.push(x);
        }
        last_read = Some(x2);
    });

    if (tag.maximum_subsequences.is_some() && tag.maximum_subsequences.unwrap() > current_sorting_bin.buffer.len()) || tag.maximum_subsequences.is_none() {
        if current_sorting_bin.buffer.len() > 10000 {
            warn!("We have a large buffer of {} reads for tag {} at iteration {}, this could take a long while...", current_sorting_bin.buffer.len(), tag.symbol, iteration);
        }
        let (buffer, correction) = current_sorting_bin.close();
        if correction.len() > largest_group_size {
            largest_group_size = correction.len();
        }

        //let shard_reader: ShardReader<SortingReadSetContainer, DefaultSort> = ShardReader::open(file).unwrap();
        for mut y in buffer {
            let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
            let corrected = correction.get(&FastaBase::to_vec_u8(&key_value.1)).unwrap();
            y.ordered_sorting_keys.push((key_value.0, FastaBase::from_vec_u8(corrected)));
            sender.send(y).unwrap();
        }
    }

    info!("For degenerate tag {} we processed {} reads, largest group size {}, we dispatched {} reads to the next level, dropped {} read collections and {} added to groups", &tag.symbol, processed_reads, largest_group_size, sent_reads, dropped_collection,added_to_group);
    sender.finished().unwrap();

    sharded_output.finish().unwrap();

    (sent_reads, ShardReader::open(aligned_temp).unwrap())
}

pub fn sort_known_level(temp_directory: &mut InstanceLivedTempDir, reader: &ShardReader<SortingReadSetContainer>, tag: &UMIConfiguration, read_count: &usize) -> (usize, ShardReader<SortingReadSetContainer>) {
    info!("Sorting known level {}",tag.symbol);
    // get the known lookup table
    info!("Loading the known lookup table for tag {}, this can take some time",tag.symbol);
    let mut known_lookup = KnownList::read_known_list_file(&tag, &tag.file.as_ref().unwrap(), &8);
    let mut processed_reads = 0;
    let mut dropped_reads = 0;
    let mut collided_reads = 0;
    info!("Sorting reads");
    let bar = ProgressBar::new(read_count.clone() as u64);

    // create a new output
    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    {
        let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&aligned_temp, 32,
                                                                                        256,
                                                                                        1 << 16).unwrap();
        let mut sender = sharded_output.get_sender();

        reader.iter_range(&Range::all()).unwrap().enumerate().for_each(|(i,x)| {
            bar.inc(1);
            processed_reads += 1;
            let mut sorting_read_set_container = x.unwrap();
            let sq = FastaBase::to_string(&sorting_read_set_container.aligned_read.aligned_read);
            let sqref = FastaBase::to_string(&sorting_read_set_container.aligned_read.aligned_ref);

            let sq2 = sorting_read_set_container.aligned_read.read_name.clone();
            assert_eq!(sorting_read_set_container.ordered_sorting_keys.len(), tag.order);

            let next_key = sorting_read_set_container.ordered_unsorted_keys.pop_front().unwrap();
            assert_eq!(next_key.0, tag.symbol);

            let corrected_hits = correct_to_known_list(&next_key.1, &mut known_lookup, &tag.max_distance);

            match (corrected_hits.hits.len(), corrected_hits.distance) {
                (x, _) if x < 1 => {
                    dropped_reads += 1;
                }
                (x, _) if x > 1 => {
                    collided_reads += 1;
                    dropped_reads += 1;
                }
                (x, y) if y > tag.max_distance => {
                    dropped_reads += 1;
                }
                (x, y) => {
                    sorting_read_set_container.ordered_sorting_keys.push((next_key.0, corrected_hits.hits.get(0).unwrap().clone()));
                    sender.send(sorting_read_set_container).unwrap();
                }
            };
        });

        info!("Dropped {} reads where we couldn't match a known barcode, {} collided reads, {} total reads",dropped_reads, collided_reads, processed_reads);
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }
    info!("For known tag {} we processed {} reads", &tag.symbol, processed_reads);

    (processed_reads - dropped_reads, ShardReader::open(aligned_temp).unwrap())
}

pub fn output_fasta_to_sorted_shard_reader(written_fasta_file: &Path,
                                           temp_directory: &mut InstanceLivedTempDir,
                                           read_structure: &SequenceLayoutDesign) -> (usize, ShardReader<SortingReadSetContainer>) {
    let aligned_temp = temp_directory.temp_file("first_pass_sorted.sharded");
    let mut processed_reads = 0;

    {
        let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(aligned_temp.clone(), 32,
                                                                                        256,
                                                                                        1 << 16).unwrap();

        let mut sender = sharded_output.get_sender();

        let input_fasta = Reader::from_file(&written_fasta_file).unwrap();
        info!("Reading fasta file {}",written_fasta_file.to_str().unwrap());
        let mut sorting_order = read_structure.umi_configurations.iter().map(|x| x.1.clone()).collect::<Vec<UMIConfiguration>>();
        sorting_order.sort_by(|a, b| a.order.cmp(&b.order));
        let sorted_tags = sorting_order.iter().map(|x| x.symbol).collect::<Vec<char>>();

        input_fasta.records().chunks(2).into_iter().for_each(|mut chunk| {
            let ref_record = chunk.nth(0).unwrap().unwrap();
            let read_record = chunk.nth(0).unwrap().unwrap();
            processed_reads += 1;
            //warn!("Reading fasta rec {}",ref_record.id().to_string());
            //warn!("Reading fasta read {}",read_record.id().to_string());

            let read_tags = extract_output_tags(&read_record.id().to_string());
            let read_tags_ordered = VecDeque::from(sorted_tags.iter().
                map(|x| (x.clone(), read_tags.1.get(x).unwrap().as_bytes().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>())).collect::<Vec<(char, Vec<FastaBase>)>>());

            let new_sorted_read_container = SortingReadSetContainer {
                ordered_sorting_keys: vec![],
                ordered_unsorted_keys: read_tags_ordered,

                aligned_read: SortedAlignment {
                    aligned_read: read_record.seq().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>(),
                    aligned_ref: ref_record.seq().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>(),
                    ref_name: "".to_string(),
                    read_name: read_tags.0,
                },
            };
            assert_eq!(new_sorted_read_container.ordered_unsorted_keys.len(), read_structure.umi_configurations.len());

            sender.send(new_sorted_read_container).unwrap();
        });
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }

    (processed_reads, ShardReader::open(aligned_temp).unwrap())
}

pub fn extract_output_tags(header_string: &String) -> (String, BTreeMap<char, String>) {
    //assert!(header_string.starts_with(">"), "The header string does not start with a > character");
    let tokens = header_string.split("_").collect::<Vec<&str>>();
    assert!(tokens.len() >= 2, "The header string does not contain any tags");

    let mut tags: BTreeMap<char, String> = BTreeMap::new();
    tokens.get(1).unwrap().split(";").for_each(|tag| {
        let tag_tokens = tag.split("=").last().unwrap().split(":").collect::<Vec<&str>>();
        assert!(tag_tokens.len() == 2, "The tag {} does not contain a key and value", tag);
        let key = char::from(tag_tokens.get(0).unwrap().as_bytes()[0].clone());
        let value = tag_tokens.get(1).unwrap().to_string();
        //println!("{} -> {}", key, value);
        tags.insert(key, value);
    });
    (String::from(tokens.get(0).unwrap().clone()).clone(), tags)
}

//pub fn order_output_keys(header_string: &String, )

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use super::*;
    use crate::read_strategies::read_set::ReadIterator;

    #[test]
    fn test_extract_output_tags() {
        let test_string = String::from(">VH00708:168:AACK37GM5:1:1101:40992:1000_key=$:TCTCACGAGGTGGCTG;key=%:CGGTTCCGAAGT;key=^:AGGGTCTCGGCC");
        let tags = extract_output_tags(&test_string);
        assert_eq!(tags.1.len(), 3);
    }

    #[test]
    fn test_consensus() {
        let basic_seqs: Vec<Vec<u8>> = vec![String::from("ATCG").as_bytes().to_vec(), String::from("GCTA").as_bytes().to_vec(), String::from("ATCG").as_bytes().to_vec()];
        let consensus = consensus(&basic_seqs);
        assert_eq!(consensus, String::from("ATCG").as_bytes().to_vec());
    }
}