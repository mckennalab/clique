use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::path::{PathBuf};
use crate::alignment_functions::{fast_align_reads};
use crate::InstanceLivedTempDir;
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration, UMISortType};
use crate::reference::fasta_reference::ReferenceManager;
use indicatif::ProgressBar;
use shardio::{Range, ShardReader, ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::{SortingReadSetContainer};
use rustc_hash::FxHashMap;
use crate::consensus::consensus_builders::write_consensus_reads;
use crate::umis::known_list::KnownList;
use crate::umis::sequence_clustering::{correct_to_known_list, get_connected_components, InputList, vantage_point_string_graph};

pub fn collapse(final_output: &String,
                fast_reference_lookup: &bool,
                temp_directory: &mut InstanceLivedTempDir,
                read_structure: &SequenceLayoutDesign,
                max_reference_multiplier: &f64,
                min_read_length: &usize,
                read1: &String,
                read2: &String,
                index1: &String,
                index2: &String,
                threads: &usize,
                find_inversions: &bool) {

    // load up the reference files
    let rm = ReferenceManager::from_yaml_input(read_structure, 12, 6);

    // validate that each reference has the specified capture groups
    let validated_references = rm.references.iter().
        map(|rf| {
            let reference_config = read_structure.references.get(String::from_utf8(rf.1.name.clone()).unwrap().as_str()).unwrap();
            SequenceLayoutDesign::validate_reference_sequence(&rf.1.sequence_u8,&reference_config.umi_configurations)
        }).all(|x| x == true);

    assert!(validated_references, "The reference sequences do not match the capture groups specified in the read structure file.");

    // align whatever combination of reads to the reference, collapsing down to a single sequence
    let aligned_temp = temp_directory.temp_file("aligned_reads.fasta");

    info!("Aligning reads to reference in directory {}", aligned_temp.as_path().display());

    let ret = fast_align_reads(&true,
                               read_structure,
                               &false,
                               &fast_reference_lookup,
                               &rm,
                               aligned_temp.as_path(),
                               max_reference_multiplier,
                               min_read_length,
                               read1,
                               read2,
                               index1,
                               index2,
                               &0.2,
                               threads,
                               find_inversions);

    info!("Sorting the aligned reads");

    let mut read_count = ret.0;

    info!("Sorting by read tags");

    // sort the reads by the tags
    let mut levels = 0;
    ret.1.into_iter().for_each(|(ref_name,sorted_reads)| {
        let mut sorted_input = sorted_reads;

        read_structure.get_sorted_umi_configurations(&ref_name).iter().for_each(|tag| {
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

        info!("writing consensus reads for reference {}",ref_name);
        // collapse the final reads down to a single sequence and write everything to the disk
        write_consensus_reads(&sorted_input, final_output, levels, &read_count, &rm, &40);
    });


}


fn consensus(input: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut consensus = Vec::new();

    // for each position
    for i in 0..input[0].len() {
        let mut counter = HashMap::new();

        // for each input string
        input.iter().for_each(|vector| {
            assert_eq!(vector.len(), input[0].len(), "string {} is not the same length as the first string {}",
                       String::from_utf8(vector.clone()).unwrap(),
                       String::from_utf8(input[0].clone()).unwrap());

            *counter.entry(&vector[i]).or_insert(0) += 1;
        });

        let mut max = 0;
        let mut consensus_byte = b'N';

        //println!("consensus {:?}",counter);
        for (byte, count) in counter {
            // if we're the new maximum OR we're tied for the maximum and we're an N or a gap, then we'll take the new value
            if count > max || (count == max && consensus_byte == b'N') || (count == max && consensus_byte == b'-') {
                max = count;
                consensus_byte = *byte;
            }
        }

        consensus.push(consensus_byte);
    }

    consensus
}

struct DegenerateBuffer {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
    writen_reads: usize,
    shard_writer: Option<Box<ShardWriter<SortingReadSetContainer>>>,
    shard_sender: Option<Box<ShardSender<SortingReadSetContainer>>>,
    output_file: PathBuf,
    tag: UMIConfiguration,
    hash_map: FxHashMap<String, usize>,
}

impl DegenerateBuffer {
    pub fn new(output_file: PathBuf, max_size: &usize, tag: UMIConfiguration) -> DegenerateBuffer {
        DegenerateBuffer {
            buffer: VecDeque::new(),
            max_buffer_size: *max_size,
            writen_reads: 0,
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
        }
    }

    /// Pushes a new item onto the buffer
    /// we buffer writing to prevent the costly disk writes when we have smaller sets of barcodes.
    /// When we overflow we write the whole buffer to disk and continue to write additional reads to
    /// disk until we've finished.
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

    /// Dumps the buffer to disk
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

    /// This function 'corrects' a list of barcodes using a connected components algorithm
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        let string_set = Vec::from_iter(self.hash_map.keys()).iter().map(|s| s.as_bytes().to_vec()).collect::<Vec<Vec<u8>>>();
        let collection = InputList { strings: string_set, max_dist: self.tag.max_distance.clone() };
        let graph = vantage_point_string_graph(&collection, self.hash_map.len() > 50000);

        let cc = get_connected_components(&graph);
        let mut final_correction: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        for group in cc {
            let minilist = InputList { strings: group, max_dist: self.tag.max_distance.clone() };

            // TODO enable this for a future improvement -- spliting connected components -- but a lot of work needs to go into validation here
            /*let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let is_subgroups = split_subgroup(&mut minigraph);

            match is_subgroups {
                None => {*/
            let group = minilist.strings.clone();
            let conc = consensus(&group); // we need to do this beforehand -- we can't do it in the loop below because we're consuming the list below which changes the results
            for s in minilist.strings {
                //info!("correct_list from {} -> {}", String::from_utf8(s.clone()).unwrap(), String::from_utf8(conc.clone()).unwrap());
                final_correction.insert(s.clone(), conc.clone());
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
        let ret = (std::mem::replace(&mut self.buffer, VecDeque::new()), final_correction);
        self.clean();
        ret
    }

    pub fn clean(&mut self) {
        self.buffer.clear();
        self.hash_map.clear();
    }
}


/// Sorts the reads by the degenerate tag
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

    let mut all_read_count: usize = 0;
    let mut output_reads: usize = 0;

    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&aligned_temp, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();

    let mut sender = sharded_output.get_sender();
    let bar = ProgressBar::new(read_count.clone() as u64);
    let mut last_read: Option<SortingReadSetContainer> = None;

    let maximum_reads_per_bin = if tag.maximum_subsequences.is_some() {
        tag.maximum_subsequences.clone().unwrap()
    } else {
        usize::MAX
    };

    let mut current_sorting_bin: Option<DegenerateBuffer> = None;

    reader.iter_range(&Range::all()).unwrap().for_each(|current_read| {
        all_read_count += 1;
        if all_read_count % 10000 == 0 {
            bar.set_position(all_read_count as u64);
        }
        let mut current_read = current_read.unwrap();
        let next_last_read = current_read.clone();

        match current_sorting_bin.as_mut() {
            None => {
                let mut bin = DegenerateBuffer::new(
                    temp_directory.temp_file(format!("{}.fasta", tag.order).as_str()),
                    &maximum_reads_per_bin,
                    tag.clone());
                bin.push(current_read);
                current_sorting_bin = Some(bin);
            }

            Some(mut bin) => {
                let reads_equal = last_read.as_ref().unwrap().cmp(&mut current_read) == Ordering::Equal;
                let exceeded_max_tags = bin.buffer.len() > maximum_reads_per_bin;

                match (reads_equal, exceeded_max_tags) {
                    (true, true) => {
                        // do nothing, we've overflowed the bin
                    }
                    (false, true) => {
                        // don't write the previous bin, but add the current read to the next bin
                        bin.clean();
                        bin.push(current_read);
                        warn!("dropping bin due to overflow, please check that this is intentional!!!")
                    }
                    (true, false) => {
                        // add the current read to the bin
                        bin.push(current_read);
                    }
                    (false, false) => {
                        // write the previous bin, and add the current read to the next bin
                        output_reads += close_and_write_bin(&mut sender, &mut bin);
                        bin.push(current_read);
                    }
                }
            }
        }

        last_read = Some(next_last_read);
    });

    match current_sorting_bin {
        None => {}
        Some(mut bin) => { output_reads += close_and_write_bin(&mut sender, &mut bin); }
    }

    bar.set_position(all_read_count as u64);

    info!("For degenerate tag {} (iteration {}) we processed {} reads, of which {} were passed to the next level", &tag.symbol, iteration, all_read_count, output_reads);
    sender.finished().unwrap();

    sharded_output.finish().unwrap();

    (output_reads, ShardReader::open(aligned_temp).unwrap())
}

fn close_and_write_bin(sender: &mut ShardSender<SortingReadSetContainer>, current_sorting_bin: &mut DegenerateBuffer) -> usize {
    let (buffer, correction) = current_sorting_bin.close();
    let ret = buffer.len();
    for mut y in buffer {
        let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
        let corrected = correction.get(&FastaBase::to_vec_u8(&key_value.1)).unwrap();
        y.ordered_sorting_keys.push((key_value.0, FastaBase::from_vec_u8(corrected)));
        sender.send(y).unwrap();
    }
    ret
}

pub fn sort_known_level(temp_directory: &mut InstanceLivedTempDir, reader: &ShardReader<SortingReadSetContainer>, tag: &UMIConfiguration, read_count: &usize) -> (usize, ShardReader<SortingReadSetContainer>) {
    info!("Sorting known level {}",tag.symbol);

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

        reader.iter_range(&Range::all()).unwrap().for_each(|x| {
            processed_reads += 1;
            if processed_reads % 10000 == 0 {
                bar.set_position(processed_reads as u64);
            }
            let mut sorting_read_set_container = x.unwrap();
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
                (_x, y) if y > tag.max_distance => {
                    dropped_reads += 1;
                }
                (_x, _y) => {
                    sorting_read_set_container.ordered_sorting_keys.push((next_key.0, corrected_hits.hits.get(0).unwrap().clone()));
                    sender.send(sorting_read_set_container).unwrap();
                }
            };
        });

        info!("Dropped {} reads (of which {} were collided reads), {} total reads",dropped_reads, collided_reads, processed_reads);
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }
    bar.set_position(processed_reads as u64);

    info!("For known tag {} we processed {} reads", &tag.symbol, processed_reads);

    (processed_reads - dropped_reads, ShardReader::open(aligned_temp).unwrap())
}

#[cfg(test)]
mod tests {
    use crate::read_strategies::read_disk_sorter::SortedAlignment;
    use crate::utils::read_utils::fake_reads;
    use super::*;

    #[test]
    fn test_consensus() {
        let basic_seqs: Vec<Vec<u8>> = vec![String::from("ATCG").as_bytes().to_vec(), String::from("GCTA").as_bytes().to_vec(), String::from("ATCG").as_bytes().to_vec()];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("ATCG").as_bytes().to_vec());

        let basic_seqs: Vec<Vec<u8>> = vec![String::from("ATCG").as_bytes().to_vec(), String::from("ATC-").as_bytes().to_vec()];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("ATCG").as_bytes().to_vec());

        // reverse the above to check that order doesn't matter (it did at one point)
        let basic_seqs: Vec<Vec<u8>> = vec![String::from("ATC-").as_bytes().to_vec(), String::from("ATCG").as_bytes().to_vec()];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("ATCG").as_bytes().to_vec());

        // real world issue
        let basic_seqs: Vec<Vec<u8>> = vec![String::from("TGGTATGCTGG-").as_bytes().to_vec(), String::from("TGGTATGCTGGG").as_bytes().to_vec()];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("TGGTATGCTGGG").as_bytes().to_vec());

        // reverse the above to check that order doesn't matter (it did at one point)
        let basic_seqs: Vec<Vec<u8>> = vec![String::from("TGGTATGCTGG-").as_bytes().to_vec(), String::from("TGGTATGCTGGG").as_bytes().to_vec()];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("TGGTATGCTGGG").as_bytes().to_vec());
    }

    #[test]
    fn test_consensus_real_world() {
        let reads = fake_reads(10, 1);
        let read_seq = reads.get(0).unwrap().read_one.seq().iter().map(|x| FastaBase::from(x.clone())).collect::<Vec<FastaBase>>();
        let fake_read = SortedAlignment {
            aligned_read: read_seq.clone(),
            aligned_ref: read_seq.clone(),
            ref_name: "".to_string(),
            read_name: "".to_string(),
            cigar_string: vec![],
            score: 0.0,
        };

        let mut tbb = DegenerateBuffer::new(
            PathBuf::from("test_data/consensus_test.fastq"),
            &1000,
            UMIConfiguration {
                symbol: 'a',
                file: None,
                reverse_complement_sequences: None,
                sort_type: UMISortType::DegenerateTag,
                length: 0,
                order: 0,
                pad: None,
                max_distance: 1,
                maximum_subsequences: None,
            });


        // real example we hit
        let st1 = SortingReadSetContainer {
            ordered_sorting_keys: vec![('a', FastaBase::from_str("AAACCCATCAGCATTA")),
                                       ('a', FastaBase::from_str("TATTGACAACCT"))],
            ordered_unsorted_keys: VecDeque::from(vec![('a', FastaBase::from_str("TGGTATGCTGG-"))]),
            aligned_read: fake_read.clone(),
        };

        let st2 = SortingReadSetContainer {
            ordered_sorting_keys: vec![('a', FastaBase::from_str("AAACCCATCAGCATTA")),
                                       ('a', FastaBase::from_str("TATTGACAACCT"))],
            ordered_unsorted_keys: VecDeque::from(vec![('a', FastaBase::from_str("TGGTATGCTGGG"))]),
            aligned_read: fake_read.clone(),
        };


        assert_eq!(st1.cmp(&st2) == Ordering::Equal, true);

        tbb.push(st1);
        tbb.push(st2.clone());
        tbb.push(st2.clone());

        let (_buffer, correction) = tbb.close();

        correction.iter().for_each(|x| {
            assert_eq!(String::from_utf8(x.1.clone()).unwrap(), String::from("TGGTATGCTGGG"));
        });
    }
}