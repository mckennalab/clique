use std::collections::{HashMap, HashSet, VecDeque};

use std::hash::BuildHasherDefault;
use std::path::PathBuf;

use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;
use collapse::LookupCollection;
use read_strategies::sequence_layout::UMISortType;
use rust_star::{DistanceGraphNode, Link, LinkedDistances, Trie};
use rustc_hash::{FxHashMap, FxHasher};
use shardio::{Range, ShardReader, ShardSender, ShardWriter};
use std::ops::Deref;
use read_strategies::read_disk_sorter::CorrectedKey;
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
    processed_sequences: usize,
}

impl SequenceCorrector {
    pub fn new(output_file: PathBuf, max_size: &usize, tag: UMIConfiguration) -> SequenceCorrector {
        SequenceCorrector {
            buffer: VecDeque::new(),
            max_buffer_size: *max_size,
            collapse_ratio: *tag.minimum_collapsing_difference.as_ref().unwrap_or(&5.0),
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
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
        if key_value.0 != self.tag.symbol {
            println!(
                "Failed read: {}\n{}\n{}\n{:?}\n{:?}\n{:?}\n{} {}",
                &item.aligned_read.read_name,
                u8s(&item.aligned_read.read_aligned),
                u8s(&item.aligned_read.reference_aligned),
                &self.tag,
                key_value,
                item.ordered_unsorted_keys,
                key_value.0,
                self.tag.symbol
            );
            println!(
                "Failed read: {}\n{}\n{}",
                &item.aligned_read.read_name,
                u8s(&item.aligned_read.read_aligned),
                u8s(&item.aligned_read.reference_aligned)
            );
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

    fn correct_known_list(&self, trie: &mut Trie) -> FxHashMap<Vec<u8>, Vec<u8>> {
        let mut knowns: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();
        let mut unmatched = 0;
        let mut multimatched = 0;
        let mut nostart = 0;
        let mut matched = 0;

        let mut sorted_tags: Vec<(Vec<u8>, usize)> = self
            .hash_map
            .iter()
            .map(|(x, y)| (x.clone(), *y))
            .collect::<Vec<(Vec<u8>, usize)>>();

        sorted_tags.sort_by(|a, b| a.0.cmp(&b.0));

        let mut search_nodes = HashSet::default();

        (0..sorted_tags.len()).for_each(|sorted_tag_index| {
            let start = if sorted_tag_index >= 1 {
                LinkedDistances::prefix_overlap_str(
                    &sorted_tags[sorted_tag_index].0,
                    &sorted_tags[sorted_tag_index - 1].0,
                )
            } else {
                0
            };
            let future = if sorted_tag_index < sorted_tags.len() - 1 {
                LinkedDistances::prefix_overlap_str(
                    &sorted_tags[sorted_tag_index + 1].0,
                    &sorted_tags[sorted_tag_index].0,
                )
            } else {
                0
            };

            if search_nodes.len() == 0 {
                search_nodes = trie.depth_links(&1);
            }
            let mut corrected_value: Vec<u8> = sorted_tags[sorted_tag_index]
                .0
                .clone()
                .into_iter()
                .filter(|x| *x != b'-')
                .collect::<Vec<u8>>();
            corrected_value.resize(self.tag.length, b'-');
            let corrected_key = corrected_value.clone();

            if start < sorted_tags[0].0.len() {
                let rt = trie.chained_search(
                    start,
                    Some(future),
                    &sorted_tags[sorted_tag_index].0,
                    &self.tag.max_distance,
                    &search_nodes,
                );
                search_nodes = rt.1;

                match rt.0.len() {
                    1 => {
                        matched += 1;
                        nostart += sorted_tags[sorted_tag_index].1;
                        knowns.insert(corrected_value, rt.0[0].0.clone());
                        debug!(
                            "{}\ttrue\t1\thit\t{}\t{}",
                            u8s(&corrected_key),
                            u8s(&rt.0[0].0),
                            sorted_tags[sorted_tag_index].1
                        );
                    }
                    0 => {
                        unmatched += 1;
                        //knowns.insert(corrected_value, sorted_tags[x].0.clone()); // be yourself dude
                        debug!(
                            "{}\tfalse\t0\tzero\tNA\t{}",
                            u8s(&corrected_key),
                            sorted_tags[sorted_tag_index].1
                        );
                    }
                    match_count => {
                        multimatched += 1;
                        // is there a minimal hit?
                        let mut min = usize::MAX;
                        let mut min_index = 0;
                        let mut min_count = 0;
                        rt.0.iter().enumerate().for_each(|(index, (_x, y))| {
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
                            nostart += sorted_tags[sorted_tag_index].1;
                            knowns.insert(corrected_value, rt.0[min_index].0.clone());
                            debug!(
                                "{}\ttrue\t{}\tmulti\tNA\t{}",
                                u8s(&corrected_key),
                                match_count,
                                sorted_tags[sorted_tag_index].1
                            );
                        } else {
                            debug!(
                                "{}\tfalse\t{}\tmulti\tNA\t{}",
                                u8s(&corrected_key),
                                match_count,
                                sorted_tags[sorted_tag_index].1
                            );
                        }
                    }
                }
            } else {
            }
        });
        debug!(
            "matched {} Unmatched {} multimatched {} nostart {} total {} super total {}",
            matched,
            unmatched,
            multimatched,
            nostart,
            self.hash_map.len(),
            self.processed_sequences
        );

        knowns
    }

    /// This function 'corrects' a list of barcodes using our Starcode clone
    pub fn correct_degenerate_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        match self.hash_map.len() {
            0 => {
                // case 0 -- no records, do nothing
                FxHashMap::default()
            }
            1 => {
                let mut knowns: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

                // TODO: wrong for known list
                // case 1 -- manually create the known list -- pad if too short
                let mut kn = self.hash_map.iter().next().unwrap().0.clone();
                if kn.len() < self.tag.length {
                    debug!("resize {} {}", self.tag.length, u8s(&kn));
                    kn.resize(self.tag.length, b'-');
                    debug!("resized {} {}", self.tag.length, u8s(&kn));
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


                let result: HashMap<Vec<u8>,Link<DistanceGraphNode>> = LinkedDistances::cluster_string_vector_list(
                    &max_length,
                    tags,
                    &self.tag.max_distance,
                    &self.collapse_ratio,
                ).into_iter().collect();

                let valid_clusters: Vec<_> = result.iter().filter(|x| x.1.borrow().valid).collect();

                let mut knowns: FxHashMap<Vec<u8>, Vec<u8>> = HashMap::default();
                let mut barcodes_to_resolve : VecDeque<Vec<u8>> = VecDeque::new();

                valid_clusters.into_iter().for_each(|(_center, dist_graph)| {
                    let connected_node = dist_graph.borrow_mut();
                    let connected_node: &DistanceGraphNode = connected_node.deref().to_owned();
                    let string_name = connected_node.string.clone();
                    knowns.insert(string_name.clone(), string_name.clone());
                    connected_node.swallowed_links.iter().for_each(|(x, y)| {
                        barcodes_to_resolve.push_front(x.clone());
                        knowns.insert(x.clone(), string_name.clone());
                    });

                    while barcodes_to_resolve.len() > 0 {
                        let processing_code = barcodes_to_resolve.pop_front().unwrap();
                        let connected_node = result.get(&processing_code).unwrap().borrow_mut();
                        let connected_node: &DistanceGraphNode = connected_node.deref().to_owned();

                        connected_node.swallowed_links.iter().for_each(|(x, y)| {
                            barcodes_to_resolve.push_front(x.clone());
                            knowns.insert(x.clone(), string_name.clone());
                        });
                    }
                });
                knowns
            }
        }
    }

    pub fn add_corrected(
        &self,
        final_correction: &FxHashMap<Vec<u8>, Vec<u8>>,
        mut y: SortingReadSetContainer,
        sender: &mut ShardSender<SortingReadSetContainer>,
    ) {
        let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
        let mut key_to_be_corrected: Vec<u8> = key_value
            .1
            .clone()
            .into_iter()
            .filter(|x| *x != b'-')
            .collect::<Vec<u8>>();
        key_to_be_corrected.resize(self.tag.length, b'-');

        let corrected = match (self.tag.sort_type, final_correction.get(&key_to_be_corrected)) {
            (UMISortType::DegenerateTag, None) => {
                panic!(
                    "Unable to find match for key {} in corrected values ({}, {}, {})",
                    u8s(&key_to_be_corrected),
                    u8s(&key_value.1),
                    final_correction.get(&key_to_be_corrected).is_some(),
                    final_correction.get(&key_value.1).is_some(),
                );
            }
            (UMISortType::KnownTag, None) => None,
            (_, Some(x)) => Some(x.clone()),
        };
        match corrected {
            Some(cv) => {
                y.ordered_sorting_keys.push((key_value.0, CorrectedKey{
                    key: self.tag.symbol,
                    original: key_to_be_corrected,
                    corrected: cv,
                }));
                sender.send(y).unwrap();
            }
            None => {}
        }
    }

    pub fn close_degenerate_list(
        &mut self,
        sender: &mut ShardSender<SortingReadSetContainer>,
    ) -> usize {
        let final_correction = self.correct_degenerate_list();
        self.close_and_write_to_shard_writer(sender, final_correction)
    }

    pub fn close_trie_known_list(
        &mut self,
        sender: &mut ShardSender<SortingReadSetContainer>,
        tag: &UMIConfiguration,
        lookup_collection: &mut LookupCollection,
    ) -> usize {
        
        match tag.levenshtein_distance {
            Some(true) => {
                let mut trie = lookup_collection
                    .ret_trie
                    .get_mut(&tag.file.clone().unwrap().clone());
                let trie = trie
                    .as_mut()
                    .unwrap();
                let final_correction = self.correct_known_list(trie);
                self.close_and_write_to_shard_writer(sender, final_correction)
            }
            None | Some(false) => {
                panic!("Calling a trie when you should of called a known list")
            }
        }
        
    }

    pub fn close_hamming_known_list(
        &mut self,
        sender: &mut ShardSender<SortingReadSetContainer>,
        tag: &UMIConfiguration,
        lookup_collection: &mut LookupCollection,
    ) -> usize {
        match tag.levenshtein_distance {
            Some(false) => {
                let mut kl = lookup_collection
                    .ret_known_lookup
                    .get_mut(&tag.file.clone().unwrap().clone());
                let kl = kl
                    .as_mut()
                    .unwrap();
                let final_correction = kl.correct_all(
                    &(self.hash_map.iter().map(|(ky,_vl)| ky.clone()).collect::<Vec<Vec<u8>>>()),
                    &(self.tag.max_distance as u32),
                );
                info!("Closing and writing corrections");
                self.close_and_write_to_shard_writer(sender, final_correction)
            }
            None | Some(true) => {
                panic!("Calling a trie when you should of called a known list")
            }
        }
        
    }

    fn close_and_write_to_shard_writer(
        &mut self,
        sender: &mut ShardSender<SortingReadSetContainer>,
        final_correction: FxHashMap<Vec<u8>, Vec<u8>>,
    ) -> usize {
        let mut read_count: usize = 0;
        let mut unbuffered_reads = 0;

        self.buffer.iter().for_each(|y| {
            self.add_corrected(&final_correction, y.clone(), sender);
            read_count += 1;
            unbuffered_reads += 1;
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
                    let current_read: SortingReadSetContainer = current_read.unwrap();
                    self.add_corrected(&final_correction, current_read, sender);

                    read_count += 1;
                    unbuffered_reads += 1;
                });
        }

        // clear everything out
        self.shard_sender = None;
        self.shard_writer = None;
        self.buffer.clear();
        self.hash_map.clear();
        debug!(
            "COUNTS {} {} {}",
            read_count,
            unbuffered_reads,
            self.hash_map.len()
        );
        read_count
    }
}

#[allow(dead_code)]
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
    use read_strategies::sequence_layout::UMISortType;
    use tempfile::NamedTempFile;

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
            levenshtein_distance: None,
        };

        // above the threshold for merging
        let tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&10, &path, &config);
        let list = tag_buffer.correct_degenerate_list();
        
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
        let tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&3, &path, &config);
        let list = tag_buffer.correct_degenerate_list();

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
        let list = tag_buffer.correct_degenerate_list();
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
        let mut tag_buffer = SequenceCorrector::new(path.clone(), &5000, config.clone());

        for _ in 0..*count {
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

        for _ in 0..*count {
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
