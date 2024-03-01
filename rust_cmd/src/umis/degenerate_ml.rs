use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;
use rustc_hash::FxHashMap;
use shardio::{ShardSender, ShardWriter};
use vpsearch::Tree;
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;
use crate::umis::known_list::FastaString;
use crate::umis::sequence_clustering::{get_connected_components, InputList, vantage_point_string_graph};

pub struct DegenerateBuffer {
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

    pub fn bin_size_exceeded(&self) -> bool {
        self.buffer.len() >= self.max_buffer_size
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
        *self
            .hash_map
            .entry(FastaBase::string(&key_value.1))
            .or_insert(0) += 1;

        match (
            &self.shard_writer,
            self.writen_reads >= self.max_buffer_size,
        ) {
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
        self.shard_writer = Some(Box::new(
            ShardWriter::new(&self.output_file, 32, 256, 1 << 16).unwrap(),
        ));
        self.shard_sender = Some(Box::new(self.shard_writer.as_mut().unwrap().get_sender()));
        self.buffer.iter().for_each(|x| {
            self.shard_sender.as_mut().unwrap().send(x.clone()).unwrap();
        });
    }

    /// This function 'corrects' a list of barcodes using a connected components algorithm
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        let string_set = Vec::from_iter(self.hash_map.keys())
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect::<Vec<Vec<u8>>>();
        let collection = InputList {
            strings: string_set,
            max_dist: self.tag.max_distance.clone(),
        };
        let graph = vantage_point_string_graph(&collection, self.hash_map.len() > 50000);

        let cc = get_connected_components(&graph);
        let mut final_correction: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        for group in cc {
            let minilist = InputList {
                strings: group,
                max_dist: self.tag.max_distance.clone(),
            };

            // TODO enable this for a future improvement -- spliting connected components -- but a lot of work needs to go into validation here
            /*let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let is_subgroups = split_subgroup(&mut minigraph);

            match is_subgroups {
                None => {*/
            let group = minilist.strings.clone();
            let conc = crate::collapse::consensus(&group); // we need to do this beforehand -- we can't do it in the loop below because we're consuming the list below which changes the results
            for s in minilist.strings {
                //info!("correct_list from {} -> {}", String::from_utf8(s.clone()).unwrap(), String::from_utf8(conc.clone()).unwrap());
                final_correction.insert(s.clone(), conc.clone());
            }
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

    pub fn close(
        &mut self,
    ) -> (
        VecDeque<SortingReadSetContainer>,
        FxHashMap<Vec<u8>, Vec<u8>>,
    ) {
        let final_correction = self.correct_list();
        let ret = (
            std::mem::replace(&mut self.buffer, VecDeque::new()),
            final_correction,
        );
        self.clean();
        ret
    }

    pub fn clean(&mut self) {
        self.buffer.clear();
        self.hash_map.clear();
    }
}

trait ListCorrector {
    fn correct_list(&self, counts: &HashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>>;
    fn get_max_distance(&self) -> u32;
}

struct MLCorrector {
    max_lev_distance: u32,
}

impl MLCorrector {
    fn find_closest(entry: &FastaString, others: &Vec<FastaString>, max_dist: u32) -> (Vec<FastaString>, u32) {
        let mut best_hits = Vec::new();
        let mut best_dist: u32 = max_dist;
        //println!("entry {} best dist {} others size {}",FastaBase::string(&entry.fa),best_dist,others.len());
        for other in others.iter() {
            let dist = entry.hamming_distance(other);
            if dist < best_dist {
                best_hits = vec![other.clone()];
                best_dist = dist;
            } else if dist == best_dist {
                best_hits.push(other.clone());
            }
            if best_dist == 0 {
                return((best_hits, best_dist));
            }
        }
        (best_hits, best_dist)
    }
}

impl ListCorrector for MLCorrector {

    fn correct_list(&self, counts: &HashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>> {
        let keys: Vec<FastaString> = counts.iter().map(|(key, value)|
            FastaString::new(FastaBase::from_str(key.as_str()))
        ).collect();

        let mut count_vec: Vec<(&String, &usize)> = counts.iter().collect();
        count_vec.sort_by(|a, b| b.1.cmp(a.1));

        let mut current_mapping: HashMap<Vec<u8>, Vec<u8>> = HashMap::default();
        let mut incorporated_list = Vec::new();

        for (fb_string, count) in count_vec {

            let vec_u8_version = fb_string.clone().into_bytes();
            let fasta_string = FastaString::new(FastaBase::from_string(fb_string));
            let best_merge_candidates = MLCorrector::find_closest(&fasta_string, &incorporated_list, self.max_lev_distance);

            if best_merge_candidates.0.is_empty() {
                current_mapping.insert(vec_u8_version.clone(), vec_u8_version.clone());
                incorporated_list.push(FastaString::from_string(fb_string));
            } else if best_merge_candidates.0.len() == 1 {
                let best = best_merge_candidates.0.get(0).unwrap();
                current_mapping.insert(vec_u8_version.clone(), FastaBase::vec_u8(&best.fa));
            } else {
                panic!("collide! {} {} {:?}",fb_string, best_merge_candidates.0.len(), best_merge_candidates.0.iter().map(|x| FastaBase::string(&x.fa)).collect::<Vec<String>>() );
            }
        }

        current_mapping
    }


    fn get_max_distance(&self) -> u32 {
        self.max_lev_distance
    }
}


#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use itertools::Itertools;
    use rand::prelude::SliceRandom;
    use crate::umis::degenerate_ml::{ListCorrector, MLCorrector};

    fn create_fake_data(min_diff: usize, length: usize) -> Vec<Vec<u8>> {
        assert_eq!(length % min_diff, 0);

        let bases = vec![vec![b'A'; min_diff], vec![b'C'; min_diff], vec![b'G'; min_diff], vec![b'T'; min_diff]];
        let mut accumulator = bases.clone();
        let base_len = bases.len();
        for i in 1..(length / min_diff) {
            let mut future_accumulator = Vec::with_capacity(accumulator.len() * base_len);
            for acc in accumulator {
                bases.iter().for_each(|base| {
                    let mut res = acc.clone();
                    res.extend(base.clone());
                    future_accumulator.push(res);
                });
            }
            accumulator = future_accumulator;
        }

        accumulator
    }


    pub fn mutate_x_positions_vec_u8(old_vec: &String, pos: usize) -> String {
        let mut rng = rand::thread_rng();

        let new_base = match old_vec.as_bytes()[pos] {
            b'A' => {vec![b'C',b'G',b'T'].choose(&mut rng).unwrap().clone()},
            b'C' => {vec![b'A',b'G',b'T'].choose(&mut rng).unwrap().clone()},
            b'G' => {vec![b'A',b'C',b'T'].choose(&mut rng).unwrap().clone()},
            b'T' => {vec![b'A',b'C',b'G'].choose(&mut rng).unwrap().clone()},
            _ => {panic!("Wrong base for {}",old_vec)}
        };

        let mut new_str = old_vec.clone();
        new_str.replace_range(pos..pos+1, String::from_utf8(vec![new_base.clone()]).unwrap().as_str());
        new_str
    }
    pub fn one_off_errors(input: &HashMap<String,usize>) -> (HashMap<String,usize>,HashMap<String,String>) {
        let mut new_hash = HashMap::new();
        let mut new_hash_mapping = HashMap::new();
        input.iter().for_each(|(k,v)| {
            new_hash.insert(k.clone(),v.clone());
            new_hash_mapping.insert(k.clone(),k.clone());

            let mutated = mutate_x_positions_vec_u8(&k,1);
            new_hash.insert(mutated.clone(), 1);
            new_hash_mapping.insert(mutated,k.clone());
        });
        (new_hash,new_hash_mapping)
    }

    #[test]
    fn test_ml_corrector() {
        let fake_list = HashMap::from_iter(create_fake_data(3, 24).into_iter().enumerate().map(|(i, c)| (String::from_utf8(c).unwrap(), i+100)));
        let mutated_lists = one_off_errors(&fake_list);
        let ml = MLCorrector { max_lev_distance: 1 };
        println!("correct {}",mutated_lists.0.len());
        let corrected = ml.correct_list(&mutated_lists.0);
        println!("correct {}",mutated_lists.0.len());
        assert_eq!(corrected.len(), fake_list.len()*2);
        //for (og, correct) in corrected {
        //    assert_eq!(mutated_lists.1.get(&String::from_utf8(og).unwrap()).unwrap().as_bytes(), correct);
        //}
    }
}

