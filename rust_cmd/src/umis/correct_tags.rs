use std::collections::{HashMap, HashSet, VecDeque};

use std::hash::BuildHasherDefault;
use std::path::PathBuf;

use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;
use bstr::ByteSlice;
use read_strategies::sequence_layout::UMISortType;
use rust_star::{DistanceGraphNode, LinkedDistances, Trie};
use rustc_hash::{FxHashMap, FxHasher};
use shardio::{Range, ShardReader, ShardSender, ShardWriter};
use std::ops::Deref;
use utils::read_utils::{strip_gaps, u8s};

pub struct SequenceCorrector {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
    collapse_ratio: f64,
    shard_writer: Option<Box<ShardWriter<SortingReadSetContainer>>>,
    shard_sender: Option<Box<ShardSender<SortingReadSetContainer>>>,
    output_file: PathBuf,
    tag: UMIConfiguration,
    hash_map: FxHashMap<Vec<u8>, usize>,
    known_tags: Option<Vec<Vec<u8>>>,
    processed_sequences: usize,
}

impl SequenceCorrector {
    pub fn new(
        output_file: PathBuf,
        max_size: &usize,
        tag: UMIConfiguration,
        known_tags: Option<&Vec<Vec<u8>>>,
    ) -> SequenceCorrector {
        SequenceCorrector {
            buffer: VecDeque::new(),
            max_buffer_size: *max_size,
            collapse_ratio: *tag.minimum_collapsing_difference.as_ref().unwrap_or(&5.0),
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
            known_tags: known_tags.map(|x| x.clone()),
            processed_sequences: 0,
        }
    }

    /// Pushes a new item onto the buffer
    /// we buffer writing to disk, to prevent the costly disk writes when we have smaller sets of barcodes.
    /// When we overflow we write the whole buffer to disk and continue to write additional reads to
    /// disk until we've finished.
    ///
    pub fn push(&mut self, mut item: SortingReadSetContainer) {
        self.processed_sequences += 1;
        assert!(self.tag.length >= self.tag.max_distance);

        let key_value = item.ordered_unsorted_keys.pop_front().unwrap();
        if(key_value.0 != self.tag.symbol) {
            
            println!("Failed read: {}\n{}\n{}\n{:?}\n{:?}\n{:?}\n{} {}",
                     &item.aligned_read.read_name,
                     u8s(&item.aligned_read.read_aligned),
                     u8s(&item.aligned_read.reference_aligned),
                     &self.tag,key_value,item.ordered_unsorted_keys,
                     key_value.0,
                     self.tag.symbol);
            println!("Failed read: {}\n{}\n{}",&item.aligned_read.read_name,u8s(&item.aligned_read.read_aligned),u8s(&item.aligned_read.reference_aligned));
            panic!("unable to process read");
        }

        // we need to make a copy of the key value, because we're going to pop it off the front of the record when we correct it in a second pass
        item.ordered_unsorted_keys.push_front(key_value.clone());

        let gapless = strip_gaps(&key_value.1);

        if gapless.len() >= self.tag.length - self.tag.max_distance
            && gapless.len() <= self.tag.length + self.tag.max_distance
        {
            *self.hash_map.entry(gapless).or_insert(0) += 1;

            // do we have a disk writer already open? or do we need to open one?
            match (
                &self.shard_sender.is_some(),
                self.buffer.len() >= self.max_buffer_size,
            ) {
                (false, true) => {
                    // we don't have a disk writer open, but we have reached the buffer size, so we need to dump the buffer to disk
                    self.buffer.push_back(item);
                    self.dump_buffer_to_disk();
                }
                (false, false) => {
                    // we don't have a disk writer open, so we can just push the item onto the buffer
                    self.buffer.push_back(item);
                }
                (true, _) => {
                    // we have a disk writer open, so we can just send the item to the disk writer
                    self.shard_sender.as_mut().unwrap().send(item).unwrap();
                }
            }
        } else {
            debug!("Dropping record *** {}", u8s(&gapless));
        }
    }
    /// Dumps the buffer to disk
    pub fn dump_buffer_to_disk(&mut self) {
        self.shard_writer = Some(Box::new(
            ShardWriter::new(&self.output_file, 32, 256, 1 << 16).unwrap(),
        ));
        self.shard_sender = Some(Box::new(self.shard_writer.as_mut().unwrap().get_sender()));
        let sender = self.shard_sender.as_mut().unwrap();
        self.buffer.iter().for_each(|x| {
            sender.send(x.clone()).unwrap();
        });
        self.buffer.clear();
    }

    /// This function 'corrects' a list of barcodes using our Starcode clone
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        let mut knowns: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();
        let mut unmatched = 0;
        let mut multimatched = 0;
        let mut nostart = 0;
        let mut matched = 0;
        match self.hash_map.len() {
            0 => {
                // case 0 -- no records, do nothing
                knowns
            }
            1 => {
                // TODO: wrong for known list
                // case 1 -- manually create the known list -- pad if too short
                let mut kn = self.hash_map.iter().next().unwrap().0.clone();
                if kn.len() < self.tag.length {
                    kn.resize(self.tag.length, b'-');
                }
                knowns.insert(kn.clone(), kn);
                knowns
            }
            _ => {
                let mut max_length: usize = 0;
                let tags = self
                    .hash_map
                    .iter()
                    .map(|x| {
                        let mut ns: Vec<u8> =
                            x.0.clone().into_iter().filter(|x| *x != b'-').collect();
                        if ns.len() < self.tag.length {
                            ns.resize(self.tag.length, b'-');
                        }
                        if ns.len() > max_length {
                            max_length = ns.len();
                        }
                        (ns, *x.1)
                    })
                    .collect::<Vec<(Vec<u8>, usize)>>();

                match self.tag.sort_type {
                    UMISortType::KnownTag => {
                        let mut trie = Trie::new(max_length);
                        //let mut file = File::create("output.txt").unwrap(); // overwrites if file exists

                        self.known_tags.as_ref().unwrap().iter().for_each(|x| {
                            trie.insert(x.as_bytes(), None, &self.tag.max_distance);
                        });

                        let mut sorted_tags: Vec<(Vec<u8>, usize)> = self
                            .hash_map
                            .iter()
                            .map(|(x, y)| (x.clone(), *y))
                            .collect::<Vec<(Vec<u8>, usize)>>();
                        sorted_tags.sort_by(|a, b| a.0.cmp(&b.0));

                        let mut search_nodes = HashSet::default();

                        (0..sorted_tags.len()).for_each(|x| {
                            let start = if x >= 1 {
                                LinkedDistances::prefix_overlap_str(
                                    &sorted_tags[x].0,
                                    &sorted_tags[x - 1].0,
                                )
                            } else {
                                0
                            };
                            let mut future = if x < sorted_tags.len() - 1 {
                                LinkedDistances::prefix_overlap_str(
                                    &sorted_tags[x + 1].0,
                                    &sorted_tags[x].0,
                                )
                            } else {
                                0
                            };

                            if search_nodes.len() == 0 {
                                search_nodes = trie.depth_links(&1);
                            }

                            if start < sorted_tags[0].0.len() {
                                let rt = trie.chained_search(
                                    start,
                                    Some(future),
                                    &sorted_tags[x].0,
                                    &self.tag.max_distance,
                                    &search_nodes,
                                );
                                search_nodes = rt.1;
                                if future < 1 {
                                    future = 1;
                                }

                                match rt.0.len() {
                                    1 => {
                                        matched += 1;
                                        nostart += sorted_tags[x].1;
                                        knowns.insert(sorted_tags[x].0.clone(), rt.0[0].0.clone());
                                        //writeln!(file, "{}\ttrue\t1\thit\t{}\t{}",u8s(&sorted_tags[x].0),u8s(&rt.0[0].0),sorted_tags[x].1).unwrap();
                                    }
                                    0 => {
                                        unmatched += 1;
                                        //writeln!(file, "{}\tfalse\t0\tzero\tNA\t{}",u8s(&sorted_tags[x].0),sorted_tags[x].1).unwrap();
                                    }
                                    x => {
                                        multimatched += 1;
                                        // is there a minimal hit?
                                        let mut min = usize::MAX;
                                        let mut min_index = 0;
                                        let mut min_count = 0;
                                        rt.0.iter().enumerate().for_each(|(index, (x, y))| {
                                            if *y < min {
                                                min = *y;
                                                min_index = index;
                                                min_count = 1;
                                            } else if *y == min {
                                                min_count += 1;
                                            }
                                        });
                                        if min_count == 1 {
                                            matched += 1;
                                            nostart += sorted_tags[x].1;
                                            knowns.insert(
                                                sorted_tags[x].0.clone(),
                                                rt.0[min_index].0.clone(),
                                            );
                                            //writeln!(file, "{}\ttrue\t{}\tmulti\tNA\t{}",u8s(&sorted_tags[x].0),x,sorted_tags[x].1).unwrap();
                                        } else {
                                            //writeln!(file, "{}\tfalse\t{}\tmulti\tNA\t{}", u8s(&sorted_tags[x].0), x, sorted_tags[x].1).unwrap();
                                        }
                                    }
                                }
                            } else {
                            }
                        });
                        println!("matched {} Unmatched {} multimatched {} nostart {} total {} super total {}",matched,unmatched,multimatched,nostart, self.hash_map.len(), self.processed_sequences);

                        knowns
                    }
                    UMISortType::DegenerateTag => {
                        let correction = LinkedDistances::cluster_string_vector_list(
                            &max_length,
                            tags,
                            &self.tag.max_distance,
                            &self.collapse_ratio,
                        );

                        correction
                            .into_iter()
                            .for_each(|(center, dist_graph)| {
                                let connected_node = dist_graph.borrow_mut();
                                let connected_node: &DistanceGraphNode =
                                    connected_node.deref().to_owned();
                                let string_name = connected_node.string.clone();
                                knowns.insert(string_name.clone(), string_name.clone());
                                //println!("Known {} to known {}",u8s(&string_name),u8s(&string_name));
                                connected_node.swallowed_links.iter().for_each(|(x, y)| {
                                    knowns.insert(x.clone(), string_name.clone());
                                    //println!("Unknown {} to known {}",u8s(&x),u8s(&string_name));
                                });
                            });
                        knowns
                    }
                }
            }
        }
    }

    pub fn close_and_write_to_shard_writer(
        &mut self,
        sender: &mut ShardSender<SortingReadSetContainer>,
    ) -> usize {
        let mut read_count: usize = 0;
        let mut buffered_reads = 0;
        let mut unbuffered_reads = 0;

        // only output status if we've looked at 30K or more outcomes
        if self.known_tags.is_some() && self.known_tags.as_ref().unwrap().len() > 30000 {
            info!("Correcting reads...");
        }

        let final_correction = self.correct_list();
        let mut hit_count = 0;
        self.buffer.iter().for_each(|y| {
            let mut y = y.clone();
            let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
            let mut corrected_value: Vec<u8> = key_value
                .1
                .clone()
                .into_iter()
                .filter(|x| *x != b'-')
                .collect::<Vec<u8>>();
            corrected_value.resize(self.tag.length, b'-');
            let corrected = match final_correction.get(&corrected_value) {
                None => {
                    println!("Correcting reads...{:?}",final_correction);
                    println!("Failed read: {}\n{}\n{}\n{:?}\n{:?}\n{:?}",
                             &y.aligned_read.read_name,
                             u8s(&y.aligned_read.read_aligned),
                             u8s(&y.aligned_read.reference_aligned),&self.tag,key_value,y.ordered_unsorted_keys);
                    panic!(
                        "Unable to find match for key {} in corrected values ({}, {}, {})",
                        u8s(&corrected_value),
                        u8s(&key_value.1),
                        final_correction.get(&corrected_value).is_some(),
                        final_correction.get(&key_value.1).is_some(),
                    );
                }
                Some(x) => x.clone(),
            };
            y.ordered_sorting_keys.push((key_value.0, corrected));
            sender.send(y).unwrap();
            read_count += 1;
            buffered_reads += 1;
        });

        if self.shard_writer.is_some() {
            self.shard_sender.as_mut().unwrap().finished().unwrap();
            self.shard_writer
                .as_mut()
                .unwrap()
                .as_mut()
                .finish()
                .unwrap();
            let reader: ShardReader<SortingReadSetContainer> =
                ShardReader::open(&self.output_file).unwrap();

            reader
                .iter_range(&Range::all())
                .unwrap()
                .for_each(|current_read| {
                    let mut current_read: SortingReadSetContainer = current_read.unwrap();
                    let key_value = current_read.ordered_unsorted_keys.pop_front().unwrap();
                    let corrected_value: Vec<u8> = key_value
                        .1
                        .clone()
                        .into_iter()
                        .filter(|x| *x != b'-')
                        .collect::<Vec<u8>>();
                    //let corrected_value = strip_gaps(&corrected_value);
                    match final_correction.get(&strip_gaps(&corrected_value)) {
                        None => {
                            //info!("Unable to find match for key {} in corrected values {} {}", u8s(&corrected_value), final_correction.contains_key(&corrected_value), u8s(&key_value.1));
                        }
                        Some(x) => {
                            hit_count += 1;
                            current_read
                                .ordered_sorting_keys
                                .push((key_value.0, x.clone()));
                            sender.send(current_read).unwrap();
                        }
                    };

                    read_count += 1;
                    unbuffered_reads += 1;
                });
        }
        // only output status if we've looked at 30K or more outcomes
        if self.known_tags.is_some() && self.known_tags.as_ref().unwrap().len() > 30000 {
            info!("Done correcting reads...");
        }

        // clear everything out
        self.shard_sender = None;
        self.shard_writer = None;
        self.buffer.clear();
        self.hash_map.clear();
        //println!("COUNTS {} {} {} {}",read_count,buffered_reads,unbuffered_reads,hit_count);
        read_count
    }
}

trait ListCorrector {
    fn correct_list(
        &self,
        counts: &FxHashMap<String, usize>,
    ) -> HashMap<Vec<u8>, Vec<u8>, BuildHasherDefault<FxHasher>>;
    fn get_max_distance(&self) -> u32;
}

#[cfg(test)]
mod tests {
    use super::*;
    use alignment::alignment_matrix::AlignmentResult;
    use bio::io::fastq::Record;
    use read_strategies::read_set::ReadSetContainer;
    use read_strategies::sequence_layout::UMISortType;
    use tempfile::tempfile;
    use tempfile::NamedTempFile;
    use {FASTA_A, FASTA_T};

    #[test]
    fn test_tag_buffer_corrects() {
        let tempfile = NamedTempFile::new().unwrap(); // creates a new temp file
        let path: PathBuf = tempfile.path().to_path_buf(); // clone the path
                                                           //let path_str = path.to_string_lossy(); // convert to a string-like type
                                                           //
                                                           //
        let config = UMIConfiguration {
            symbol: '1',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 10,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: Some(5000),
            minimum_collapsing_difference: Some(5.0),
            max_gaps: None,
        };

        // above the threshold for merging
        let mut tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&10, &path, &config);
        let list = tag_buffer.correct_list();
        let minimum_thresh = 5.0;

        assert!(list.contains_key("AAAAATTTTT".as_bytes()));
        assert_eq!(
            list.get("AAAAATTTTT".as_bytes()).unwrap().clone(),
            "AAAAATTTTT".as_bytes().to_vec()
        );

        assert!(list.contains_key("AAAAATTTGT".as_bytes()));
        assert_eq!(
            list.get("AAAAATTTGT".as_bytes()).unwrap().clone(),
            "AAAAATTTTT".as_bytes().to_vec()
        );

        assert!(list.contains_key("GGGGGCCCCC".as_bytes()));
        assert_eq!(
            list.get("GGGGGCCCCC".as_bytes()).unwrap().clone(),
            "GGGGGCCCCC".as_bytes().to_vec()
        );

        assert!(list.contains_key("GCGGGCCCCC".as_bytes()));
        assert_eq!(
            list.get("GCGGGCCCCC".as_bytes()).unwrap().clone(),
            "GGGGGCCCCC".as_bytes().to_vec()
        );

        // below the threshold for merging
        let mut tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&3, &path, &config);
        let list = tag_buffer.correct_list();

        assert!(list.contains_key("AAAAATTTTT".as_bytes()));
        assert_eq!(
            list.get("AAAAATTTTT".as_bytes()).unwrap().clone(),
            "AAAAATTTTT".as_bytes().to_vec()
        );

        assert!(list.contains_key("AAAAATTTGT".as_bytes()));
        assert_eq!(
            list.get("AAAAATTTGT".as_bytes()).unwrap().clone(),
            "AAAAATTTGT".as_bytes().to_vec()
        ); // not corrected

        assert!(list.contains_key("GGGGGCCCCC".as_bytes()));
        assert_eq!(
            list.get("GGGGGCCCCC".as_bytes()).unwrap().clone(),
            "GGGGGCCCCC".as_bytes().to_vec()
        );

        assert!(list.contains_key("GCGGGCCCCC".as_bytes()));
        assert_eq!(
            list.get("GCGGGCCCCC".as_bytes()).unwrap().clone(),
            "GCGGGCCCCC".as_bytes().to_vec()
        ); // not corrected

        let mut tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&10, &path, &config);
        let off_by_one_1 = "GGGGGCCCC-";
        let off_by_two_1 = "GGGGGCCCCA";
        let off_by_two_2 = "GGGGCCCCC-";
        let truth_seqeun = "GGGGGCCCCC";

        tag_buffer.push(create_fake_read_set_container(
            &"read1".to_string(),
            &off_by_one_1.to_string(),
            &config,
        ));
        tag_buffer.push(create_fake_read_set_container(
            &"read1".to_string(),
            &off_by_two_1.to_string(),
            &config,
        ));
        tag_buffer.push(create_fake_read_set_container(
            &"read1".to_string(),
            &off_by_two_2.to_string(),
            &config,
        ));
        let list = tag_buffer.correct_list();
        list.iter()
            .for_each(|(x, y)| println!("x {} y {}", u8s(x), u8s(y)));

        assert!(list.contains_key(off_by_one_1.as_bytes()));
        assert!(list.contains_key(off_by_two_1.as_bytes()));
        assert!(list.contains_key(off_by_two_2.as_bytes()));
        assert_eq!(
            list.get(off_by_one_1.as_bytes()).unwrap().clone(),
            truth_seqeun.as_bytes().to_vec()
        ); // not corrected
        assert_eq!(
            list.get(off_by_two_1.as_bytes()).unwrap().clone(),
            truth_seqeun.as_bytes().to_vec()
        ); // not corrected
        assert_eq!(
            list.get(off_by_two_2.as_bytes()).unwrap().clone(),
            truth_seqeun.as_bytes().to_vec()
        ); // not corrected
    }

    fn create_tag_buffer_with_set_anchor_seq_count(
        count: &usize,
        path: &PathBuf,
        config: &UMIConfiguration,
    ) -> SequenceCorrector {
        let mut tag_buffer = SequenceCorrector::new(path.clone(), &5000, config.clone(), None);

        for i in 0..*count {
            tag_buffer.push(create_fake_read_set_container(
                &"read1".to_string(),
                &"AAAAATTTTT".to_string(),
                &config,
            ));
        }
        tag_buffer.push(create_fake_read_set_container(
            &"read1".to_string(),
            &"AAAAATTTGT".to_string(),
            &config,
        ));

        for i in 0..*count {
            tag_buffer.push(create_fake_read_set_container(
                &"read1".to_string(),
                &"GGGGGCCCCC".to_string(),
                &config,
            ));
        }
        tag_buffer.push(create_fake_read_set_container(
            &"read1".to_string(),
            &"GCGGGCCCCC".to_string(),
            &config,
        ));

        tag_buffer
    }

    fn create_fake_read_set_container(
        name: &String,
        read_one_seq: &String,
        config: &UMIConfiguration,
    ) -> SortingReadSetContainer {
        let read_1_qual = vec![b'H'; read_one_seq.len()];

        let mut unsorted_keys = VecDeque::default();
        unsorted_keys.push_front((
            config.symbol.clone(),
            read_one_seq.clone().into_bytes().to_vec(),
        ));

        SortingReadSetContainer {
            ordered_sorting_keys: Vec::new(),
            ordered_unsorted_keys: unsorted_keys,
            aligned_read: AlignmentResult {
                reference_name: "default".to_string(),
                read_name: name.clone(),
                reference_aligned: vec![], // can be empty -- only the tags are used in this case
                read_aligned: vec![],      // can be empty -- only the tags are used in this case
                read_quals: None,          // can be empty -- only the tags are used in this case
                cigar_string: vec![],      // can be empty -- only the tags are used in this case
                path: vec![],              // can be empty -- only the tags are used in this case
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        }
    }
}
