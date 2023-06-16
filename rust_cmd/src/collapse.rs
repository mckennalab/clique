use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::path::{Path, PathBuf};
use crate::alignment_functions::align_reads;
use crate::InstanceLivedTempDir;
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration, UMISortType};
use crate::reference::fasta_reference::ReferenceManager;
use bio::io::fasta::*;
use itertools::Itertools;
use shardio::{DefaultSort, Range, ShardReader, ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};
use crate::sequence_lookup::KnownLookup;
use rustc_hash::FxHashMap;
use crate::consensus::consensus_builders::write_consensus_reads;
use crate::umis::sequence_clustering::{get_connected_components, input_list_to_graph, InputList, split_subgroup, string_distance};

pub fn collapse(reference: &String,
                final_output: &String,
                temp_directory: &mut InstanceLivedTempDir,
                read_structure: &SequenceLayoutDesign,
                max_reference_multiplier: &usize,
                min_read_length: &usize,
                read1: &String,
                read2: &String,
                index1: &String,
                index2: &String,
                threads: &usize,
                inversions: &bool) {


    // load up the reference files
    let rm = ReferenceManager::from(&reference, 8);

    // validate that each reference has the specified capture groups
    let validated_references = rm.references.iter().
        map(|rf| read_structure.validate_reference_sequence(&rf.sequence_u8)).all(|x| x == true);

    assert!(validated_references, "The reference sequences do not match the capture groups specified in the read structure file.");

    // align whatever combination of reads to the reference, collapsing down to a single sequence
    let aligned_temp = temp_directory.temp_file("aligned_reads.fasta");

    warn!("Aligning reads to reference in directory {}", aligned_temp.as_path().display());

    align_reads(&true,
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
                threads,
                inversions);

    warn!("Sorting the aligned reads");

    // TODO: fix alignment to output as a sorted file and get rid of this rewrite step
    let mut sorted_input = output_fasta_to_sorted_shard_reader(&aligned_temp.as_path(), temp_directory, read_structure);

    warn!("Sorting by read tags");

    // sort the reads by the tags
    read_structure.get_sorted_umi_configurations().iter().for_each(|tag| {
        match tag.sort_type {
            UMISortType::KnownTag => {
                sorted_input = sort_known_level(temp_directory, &sorted_input, &tag);
            }
            UMISortType::DegenerateTag => {
                sorted_input = sort_degenerate_level(temp_directory, &sorted_input, &tag);
            }
        }
    });

    // collapse the final reads down to a single sequence and write everything to the disk
    write_consensus_reads(&sorted_input, final_output, 1);
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
        let graph = input_list_to_graph(&collection, string_distance, false);

        let cc = get_connected_components(&graph);
        println!("CC SIZE: {}", &cc.len());

        let mut final_correction: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        for group in cc {
            let minilist = InputList { strings: group, max_dist: u64::try_from(self.tag.max_distance.clone()).unwrap() };
            let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let is_subgroups = split_subgroup(&mut minigraph);

            match is_subgroups {
                None => {
                    let group = minilist.strings.clone();
                    for s in minilist.strings {
                        final_correction.insert(s.clone(), consensus(&group));
                    }
                }
                Some(x) => {
                    for sgroup in &x {
                        for s in sgroup {
                            final_correction.insert(s.clone(), consensus(&sgroup));
                        }
                    }
                }
            }
        }
        final_correction
    }

    pub fn close(&mut self) -> (PathBuf, FxHashMap<Vec<u8>, Vec<u8>>) {
        let final_correction = self.correct_list();
        let fl = self.output_file.clone();
        if self.shard_writer.is_none() {
            self.dump_buffer();
        }
        self.shard_sender.as_mut().unwrap().finished().unwrap();
        self.shard_writer.as_mut().unwrap().finish().unwrap();

        (fl, final_correction)
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
pub fn sort_degenerate_level(temp_directory: &mut InstanceLivedTempDir, reader: &ShardReader<SortingReadSetContainer>, tag: &UMIConfiguration) -> ShardReader<SortingReadSetContainer> {
    warn!("Sorting degenerate level {}",tag.symbol);
    // create a new output
    let mut processed_reads = 0;
    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&aligned_temp, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();

    let mut sender = sharded_output.get_sender();

    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut current_sorting_bin: TagStrippingDiskBackedBuffer =
        TagStrippingDiskBackedBuffer::new(temp_directory.temp_file(format!("{}_{}.fasta", 0, tag.order).as_str()), 10000, tag.clone());

    reader.iter_range(&Range::all()).unwrap().enumerate().for_each(|(i, x)| {
        let mut x = x.unwrap().clone();
        processed_reads += 1;
        if processed_reads % 10000 == 0 {
            warn!("Processed {} reads", processed_reads);
        }
        if !last_read.is_none() {
            if last_read.as_ref().unwrap().cmp(&&mut x) == Ordering::Equal {
                current_sorting_bin.push(x);
            } else {
                let (file, correction) = current_sorting_bin.close();
                let shard_reader: ShardReader<SortingReadSetContainer, DefaultSort> = ShardReader::open(file).unwrap();
                for y in shard_reader.iter_range(&Range::all()).unwrap() {
                    let mut y: SortingReadSetContainer = y.unwrap().clone();
                    let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
                    let corrected = correction.get(&FastaBase::to_vec_u8(&key_value.1)).unwrap();
                    y.ordered_sorting_keys.push((key_value.0, FastaBase::from_vec_u8(corrected)));
                    sender.send(y).unwrap();
                }

                current_sorting_bin = TagStrippingDiskBackedBuffer::new(temp_directory.temp_file(format!("{}_{}.fasta", i, tag.order).as_str()), 10000, tag.clone());
                current_sorting_bin.push(x);
            }
        } else {
            current_sorting_bin.push(x);
        }
    });
    let (file, correction) = current_sorting_bin.close();
    let shard_reader: ShardReader<SortingReadSetContainer, DefaultSort> = ShardReader::open(file).unwrap();
    for y in shard_reader.iter_range(&Range::all()).unwrap() {
        let mut y: SortingReadSetContainer = y.unwrap().clone();
        let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
        let corrected = correction.get(&FastaBase::to_vec_u8(&key_value.1)).unwrap();
        y.ordered_sorting_keys.push((key_value.0, FastaBase::from_vec_u8(corrected)));
        sender.send(y).unwrap();
    }

    warn!("For degenerate tag {} we processed {} reads", &tag.symbol, processed_reads);
    sender.finished().unwrap();

    sharded_output.finish().unwrap();

    ShardReader::open(aligned_temp).unwrap()
}

pub fn sort_known_level(temp_directory: &mut InstanceLivedTempDir, reader: &ShardReader<SortingReadSetContainer>, tag: &UMIConfiguration) -> ShardReader<SortingReadSetContainer> {
    warn!("Sorting known level {}",tag.symbol);
    // get the known lookup table
    warn!("Loading the known lookup table for tag {}, this can take some time",tag.symbol);
    let known_lookup = KnownLookup::from(tag);
    let mut processed_reads = 0;
    warn!("Sorting reads");

    // create a new output
    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    {
        let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&aligned_temp, 32,
                                                                                        256,
                                                                                        1 << 16).unwrap();
        let mut sender = sharded_output.get_sender();
        let mut dropped_reads = 0;
        reader.iter_range(&Range::all()).unwrap().for_each(|x| {
            processed_reads += 1;
            let mut sorting_read_set_container = x.unwrap();
            assert_eq!(sorting_read_set_container.ordered_sorting_keys.len(), tag.order);

            let next_key = sorting_read_set_container.ordered_unsorted_keys.pop_front().unwrap();
            assert_eq!(next_key.0, tag.symbol);

            match known_lookup.correct(&FastaBase::to_string(&next_key.1), &tag.max_distance, false) {
                None => { dropped_reads += 1 }
                Some(x) => {
                    sorting_read_set_container.ordered_sorting_keys.push((next_key.0, FastaBase::from_string(&x)));
                    assert_eq!(sorting_read_set_container.ordered_sorting_keys.len(), &tag.order + 1);
                    sender.send(sorting_read_set_container).unwrap();
                }
            }

            if processed_reads % 10000 == 0 {
                warn!("Processed {} reads", processed_reads);
            }
        });

        warn!("Dropped {} where we coudn't match a known barcode",dropped_reads);
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }
    warn!("For known tag {} we processed {} reads", &tag.symbol, processed_reads);

    ShardReader::open(aligned_temp).unwrap()
}

pub fn output_fasta_to_sorted_shard_reader(written_fasta_file: &Path,
                                           temp_directory: &mut InstanceLivedTempDir,
                                           read_structure: &SequenceLayoutDesign) -> ShardReader<SortingReadSetContainer> {
    let aligned_temp = temp_directory.temp_file("first_pass_sorted.sharded");
    {
        let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(aligned_temp.clone(), 32,
                                                                                        256,
                                                                                        1 << 16).unwrap();

        let mut sender = sharded_output.get_sender();

        let input_fasta = Reader::from_file(&written_fasta_file).unwrap();
        warn!("Reading fasta file {}",written_fasta_file.to_str().unwrap());
        let mut sorting_order = read_structure.umi_configurations.iter().map(|x| x.1.clone()).collect::<Vec<UMIConfiguration>>();
        sorting_order.sort_by(|a, b| a.order.cmp(&b.order));
        let sorted_tags = sorting_order.iter().map(|x| x.symbol).collect::<Vec<char>>();

        input_fasta.records().chunks(2).into_iter().for_each(|mut chunk| {
            let ref_record = chunk.nth(0).unwrap().unwrap();
            let read_record = chunk.nth(0).unwrap().unwrap();
            //warn!("Reading fasta rec {}",ref_record.id().to_string());
            //warn!("Reading fasta read {}",read_record.id().to_string());

            let read_tags = extract_output_tags(&read_record.id().to_string());
            let read_tags_ordered = VecDeque::from(sorted_tags.iter().
                map(|x| (x.clone(), read_tags.get(x).unwrap().as_bytes().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>())).collect::<Vec<(char, Vec<FastaBase>)>>());

            let new_sorted_read_container = SortingReadSetContainer {
                ordered_sorting_keys: vec![],
                ordered_unsorted_keys: read_tags_ordered,

                aligned_read: SortedAlignment {
                    aligned_read: read_record.seq().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>(),
                    aligned_ref: ref_record.seq().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>(),
                    ref_name: "".to_string(),
                },
            };
            assert_eq!(new_sorted_read_container.ordered_unsorted_keys.len(), read_structure.umi_configurations.len());

            sender.send(new_sorted_read_container).unwrap();
        });
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }

    ShardReader::open(aligned_temp).unwrap()
}

pub fn extract_output_tags(header_string: &String) -> BTreeMap<char, String> {
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
    tags
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
        assert_eq!(tags.len(), 3);
    }

    #[test]
    fn test_consensus() {
        let basic_seqs: Vec<Vec<u8>> = vec![String::from("ATCG").as_bytes().to_vec(), String::from("GCTA").as_bytes().to_vec(), String::from("ATCG").as_bytes().to_vec()];
        let consensus = consensus(&basic_seqs);
        assert_eq!(consensus, String::from("ATCG").as_bytes().to_vec());
    }
}