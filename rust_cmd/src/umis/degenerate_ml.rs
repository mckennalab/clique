use std::collections::{HashMap, HashSet, VecDeque};
use std::hash::BuildHasherDefault;
use std::path::PathBuf;
use std::thread::current;
use itertools::enumerate;
use num_traits::Bounded;
use rustc_hash::{FxHasher, FxHashMap};
use shardio::{ShardSender, ShardWriter};
use vpsearch::{BestCandidate, MetricSpace, Tree};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;
use crate::umis::known_list::FastaString;
use crate::umis::sequence_clustering::{get_connected_components, InputList, RadiusBasedNeighborhood, vantage_point_string_graph};

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
    fn correct_list(&self, counts: &FxHashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>, BuildHasherDefault<FxHasher>>;
    fn get_max_distance(&self) -> u32;
}

struct MLCorrector {
    max_lev_distance: u32,
}

impl MLCorrector {
    /*fn find_closest(entry: &FastaString, vantage_tree: &Tree<FastaString>, others: &Vec<FastaString>, max_dist: u32) -> (Vec<FastaString>, u32) {
        let mut best_hits = Vec::new();
        let mut best_dist: u32 = max_dist;
        //println!("entry {} best dist {} others size {}",FastaBase::string(&entry.fa),best_dist,others.len());
        for other in others.iter() {
            let dist = entry.distance(other,&()); //entry.hamming_distance(other);
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
    }*/

}


/// Add custom search for finding the index of the N nearest points
struct CountBasedNeighborhood<Item: MetricSpace<Impl>, Impl> {
    // Max amount of items
    max_item_count: usize,
    // The max distance we have observed so far
    max_observed_distance: Item::Distance,
    // A list of indexes no longer than max_item_count sorted by distance
    distance_x_index: Vec<(Item::Distance, usize)>,
    // the maximum distance we'd consider a hit
    max_distance: Item::Distance,
    strings: Vec<FastaString>,
}

impl<Item: MetricSpace<Impl>, Impl> CountBasedNeighborhood<Item, Impl> {
    /// Helper function for creating the CountBasedNeighborhood struct.
    /// Here `item_count` is the amount of items returned, the k in knn.
    fn new(item_count: usize, max_distance: Item::Distance, strings: Vec<FastaString>) -> Self {
        CountBasedNeighborhood {
            max_item_count: item_count,
            max_observed_distance: <Item::Distance as Bounded>::min_value(),
            distance_x_index: Vec::<(Item::Distance, usize)>::new(),
            max_distance,
            strings
        }
    }

    /// Insert a single index in the correct possition given that the
    /// `distance_x_index` is already sorted.
    fn insert_index(&mut self, index: usize, distance: Item::Distance) {
        // Add the new item at the end of the list.
        self.distance_x_index.push((distance, index));
        // We only need to sort lists with more than one entry
        if self.distance_x_index.len() > 1 {
            // Start indexing at the end of the vector. Note that len() is 1 indexed.
            let mut n = self.distance_x_index.len() - 1;
            // at n is further than n -1 we swap the two.
            // Prefrom a single insertion sort pass. If the distance of the element
            while n > 0 && self.distance_x_index[n].0 < self.distance_x_index[n - 1].0 {
                self.distance_x_index.swap(n, n - 1);
                n = n - 1;
            }
            //self.distance_x_index.truncate(self.max_item_count);
        }
        // Update the max observed distance, unwrap is safe because this function
        // inserts a point and the `max_item_count` is more then 0.
        self.max_observed_distance = self.distance_x_index.last().unwrap().0
    }
}

/// Best candidate definitions that tracks of the index all the points
/// within the radius of `distance` as specified in the `RadiusBasedNeighborhood`.
impl<Item: MetricSpace<Impl> + Clone, Impl> BestCandidate<Item, Impl>
for CountBasedNeighborhood<Item, Impl>
    where <Item as MetricSpace<Impl>>::Distance: std::fmt::Display
{
    type Output = Vec<(usize,<Item as MetricSpace<Impl>>::Distance)>;

    #[inline]
    fn consider(
        &mut self,
        _: &Item,
        distance: Item::Distance,
        candidate_index: usize,
        _: &Item::UserData,
    ) {
        //println!("item {} distance {} cand index {}",item,distance,candidate_index);
        // If the distance is lower than the max_observed distance we
        // need to add that index into the sorted_ids and update the
        // `max_observed_distance`. If the sorted_ids is already at max
        // capacity we drop the point with the max distance and find
        // out what the new `max_observed_distance` is by looking at
        // the last entry in the `distance_x_index` vector. We also
        // include the point if the `distance_x_index` is not full yet.
        if distance <= self.max_distance
        {
            self.insert_index(candidate_index, distance);
        }
    }

    #[inline]
    fn distance(&self) -> Item::Distance {
        // return distance of the Nth farthest as we have currently observed it.
        // All other points currently in the state will be closer than this.
        self.max_distance
    }

    fn result(self, _: &Item::UserData) -> Self::Output {
        self.distance_x_index.iter().map(|(ds,ind)| (*ind,ds.clone())).collect()
    }
}

impl ListCorrector for MLCorrector {


    fn correct_list(&self, counts: &FxHashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>, BuildHasherDefault<FxHasher>> {
        let keys: Vec<FastaString> = counts.iter().enumerate().map(|(index,(key, value))| {
            FastaString::new(FastaBase::from_str(key.as_str()))
        }).collect();

        let mut count_vec: Vec<(&String, &usize)> = counts.iter().collect();
        count_vec.sort_by(|a, b| b.1.cmp(a.1));

        let mut incorporated_list = HashMap::new();

        let vantage_tree = Tree::new(&keys);
        let mut current_mapping: FxHashMap<Vec<u8>, Vec<u8>> = HashMap::default();

        for (index, (fb_string, count)) in enumerate(count_vec) {
            let vec_u8_version = fb_string.clone().into_bytes();
            let fasta_string = FastaString::new(FastaBase::from_string(fb_string));
            if index == 0 {
                current_mapping.insert(vec_u8_version.clone(),vec_u8_version.clone());
                incorporated_list.insert(fasta_string,true);
            } else {
                let best_merge_candidates = vantage_tree.find_nearest_custom(
                    &fasta_string,
                    &(),
                    CountBasedNeighborhood::new(1, 1, keys.clone())
                );
                if best_merge_candidates.is_empty() {
                    current_mapping.insert(vec_u8_version.clone(), vec_u8_version.clone());
                    incorporated_list.insert(FastaString::from_string(fb_string), true);

                } else {
                    // find the best valid hit
                    let first_best : Vec<(usize, bool, u32)>= best_merge_candidates.into_iter().map(|(ht,dist)| {
                        let strin = keys.get(ht).unwrap().clone();
                        let is_in = incorporated_list.contains_key(&strin);
                        (ht,is_in, dist)
                    }).filter(|(p,x, dist)| *x && *dist <= self.max_lev_distance).collect();

                    match first_best.len() {
                        x if x == 0 => {
                            current_mapping.insert(vec_u8_version.clone(), vec_u8_version.clone());
                            incorporated_list.insert(fasta_string, true);
                        }
                        x => {
                            let best = keys.get(first_best.iter().next().unwrap().0).unwrap();

                            current_mapping.insert(vec_u8_version.clone(), FastaBase::vec_u8(&best.fa));
                        }
                    }

                }
                if index % 1000 == 0 {
                    println!("Processed {} tags", index);
                }
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
    use std::fs::File;
    use itertools::Itertools;
    use rand::prelude::SliceRandom;
    use rustc_hash::FxHashMap;
    use triple_accel::levenshtein;
    use vpsearch::Tree;
    use crate::alignment::fasta_bit_encoding::FastaBase;
    use crate::umis::degenerate_ml::{CountBasedNeighborhood, ListCorrector, MLCorrector};
    use crate::umis::known_list::FastaString;

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
            b'A' => b'C',//{vec![b'C',b'G',b'T'].choose(&mut rng).unwrap().clone()},
            b'C' => b'A',//{vec![b'A',b'G',b'T'].choose(&mut rng).unwrap().clone()},
            b'G' => b'T',//{vec![b'A',b'C',b'T'].choose(&mut rng).unwrap().clone()},
            b'T' => b'G',//{vec![b'A',b'C',b'G'].choose(&mut rng).unwrap().clone()},
            _ => {panic!("Wrong base for {}",old_vec)}
        };

        let mut new_str = old_vec.clone();
        new_str.replace_range(pos..pos+1, String::from_utf8(vec![new_base.clone()]).unwrap().as_str());
        new_str
    }

    pub fn one_off_errors(input: &FxHashMap<String,usize>) -> (FxHashMap<String,usize>,FxHashMap<String,String>) {
        let mut new_hash = FxHashMap::default();
        let mut new_hash_mapping = FxHashMap::default();
        input.iter().for_each(|(k,v)| {
            new_hash.insert(k.clone(),v.clone());
            new_hash_mapping.insert(k.clone(),k.clone());

            let mutated = mutate_x_positions_vec_u8(&k,1);
            new_hash.insert(mutated.clone(), 1);
            new_hash_mapping.insert(mutated,k.clone());
        });
        (new_hash,new_hash_mapping)
    }
    use std::io::Write;

    #[test]
    fn test_ml_corrector() {
        let fake_list = FxHashMap::from_iter(create_fake_data(1, 6).into_iter().enumerate().map(|(i, c)| (String::from_utf8(c).unwrap(), i+100)));
        let mutated_lists = one_off_errors(&fake_list);
        let ml = MLCorrector { max_lev_distance: 1 };
        let mut output_fl = File::create("test_output.pairs").unwrap();
        for (seq, correct) in mutated_lists.1.clone() {
            write!(output_fl,"{} {}\n",seq,correct);
        }
        println!("correct {}",mutated_lists.0.len());
        use std::time::Instant;
        let now = Instant::now();
        let corrected = ml.correct_list(&mutated_lists.0);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
        println!("correct {}",mutated_lists.0.len());
        assert_eq!(corrected.len(), fake_list.len());
        for (og, correct) in corrected {
            assert_eq!(mutated_lists.1.get(&String::from_utf8(og).unwrap()).unwrap().as_bytes(), correct);
        }
    }

    #[test]
    fn time_ml_corrector() {
        use std::time::Instant;
        for barcode_size in 4..16 {
            let fake_list = FxHashMap::from_iter(create_fake_data(1, barcode_size).into_iter().enumerate().map(|(i, c)| (String::from_utf8(c).unwrap(), i+100)));
            let keys: Vec<FastaString> = fake_list.iter().enumerate().map(|(index,(key, value))| {
                FastaString::new(FastaBase::from_str(key.as_str()))
            }).collect();
            let ml = MLCorrector { max_lev_distance: 1 };
            let now = Instant::now();
            let vantage_tree = Tree::new(&keys);
            let elapsed = now.elapsed();
            println!("Run size {} ({}) Elapsed: {:.2?}", barcode_size, fake_list.len(), elapsed);
        }
    }


    fn test_tree_construction() {
        let keys = vec![FastaString::from_string(&"TGTTTTTTTCCC".to_ascii_uppercase()),
                        FastaString::from_string(&"TTTTTTTTTCCC".to_ascii_uppercase())];

        let vantage_tree = Tree::new(&keys);

        let search = FastaString::from_string(&"TGTTTTTTTCCC".to_ascii_uppercase());

        let best_merge_candidates = vantage_tree.find_nearest_custom(
            &search,
            &(),
            CountBasedNeighborhood::new(1, 1, keys.clone()));

        assert_eq!(best_merge_candidates.len(),2);

        let search = FastaString::from_string(&"TTTTTTTTTCCC".to_ascii_uppercase());

        // check a couple hundred times
        for i in 0..500 {
            let best_merge_candidates = vantage_tree.find_nearest_custom(
                &search,
                &(),
                CountBasedNeighborhood::new(1, 1, keys.clone()));

            assert_eq!(best_merge_candidates.len(), 2);
        }
        let search = FastaString::from_string(&"AATTTTTTTCCC".to_ascii_uppercase());

        let best_merge_candidates = vantage_tree.find_nearest_custom(
            &search,
            &(),
            CountBasedNeighborhood::new(1, 1, keys.clone()));

        assert_eq!(best_merge_candidates.len(),0);


    }
}

