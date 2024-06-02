
use std::collections::{HashMap, HashSet};
use petgraph::algo::{tarjan_scc};
use petgraph::prelude::*;
use indicatif::ProgressBar;
use vpsearch::{BestCandidate, MetricSpace};
use crate::alignment::fasta_bit_encoding::FastaBase;






pub struct InputList {
    pub strings: Vec<Vec<u8>>,
    pub max_dist: usize, // the minimum distance to consider for creating an edge
}

pub struct StringGraph {
    pub graph: GraphMap<u32, u32, Undirected>,
    pub string_to_node: HashMap<Vec<u8>, u32>,
    pub node_to_string: HashMap<u32, Vec<u8>>,
    pub max_distance: u64,
}

pub fn string_distance_no_break(str1: &Vec<u8>, str2: &Vec<u8>, _max_dist: &usize) -> usize {
    //assert_eq!(str1.len(), str2.len());
    str1.iter().zip(str2.iter()).map(|(c1, c2)| if c1 == c2 { 0 } else { 1 }).sum()
}

#[allow(dead_code)]
pub fn string_distance_break(str1: &Vec<u8>, str2: &Vec<u8>, max_dist: &usize) -> usize {
    //assert_eq!(str1.len(), str2.len());
    let mut dist: usize = 0;
    for (c1, c2) in str1.iter().zip(str2.iter()) {
        if c1 != c2 {
            dist += 1;
            if &dist > max_dist {
                return dist;
            }
        }
    }
    dist
}

pub struct BestHits {
    pub hits: Vec<Vec<FastaBase>>,
    pub distance: usize,
}

impl Clone for BestHits {
    fn clone(&self) -> BestHits {
        BestHits { hits: self.hits.clone(), distance: self.distance }
    }
}

#[allow(dead_code)]
pub fn average_dist(strings: &Vec<Vec<u8>>,compare: fn(&Vec<u8>, &Vec<u8>) -> u64) -> f64 {
    let mut dist : f64 = 0.0;
    let mut count: usize = 0;

    for st1 in strings {
        for st2 in strings {
            dist += compare(st1,st2) as f64;
            count += 1;
        }
    }
    dist / (count as f64)
}



#[derive(Clone)]
struct U8String {
    u8str : Vec<u8>
}

impl MetricSpace for U8String {
    type UserData = ();
    type Distance = u32;

    fn distance(&self, other: &Self, _: &Self::UserData) -> Self::Distance {
        let mut dist: u32 = 0;
        for (c1, c2) in self.u8str.iter().zip(other.u8str.iter()) {
            if c1 != c2 {
                dist += 1;
            }
        }
        dist
    }
}

/// Add custom search for finding the index of multiple points in a radius
/// The index of all point with a euclidean distance strictly less than
/// `max_distance` will be returned.
pub struct RadiusBasedNeighborhood<Item: MetricSpace<Impl>, Impl> {
    max_distance: Item::Distance,
    ids: HashSet<u32>,
    distances: HashMap<u32,Item::Distance>,

}

impl<Item: MetricSpace<Impl>, Impl> RadiusBasedNeighborhood<Item, Impl> {
    /// Helper function for creating the RadiusBasedNeighborhood struct.
    /// Here `max_distance` is an exclusive upper bound to the euclidean distance.
    pub fn new(max_distance: Item::Distance) -> Self {
        RadiusBasedNeighborhood {
            max_distance,
            ids: HashSet::<u32>::new(),
            distances: HashMap::new(),
        }
    }
}

#[derive(Clone,Copy,Hash)]
struct ClosestHits {
    id: u32,
    distance: u32,
}

/// Best candidate definitions that tracks of the index all the points
/// within the radius of `distance` as specified in the `RadiusBasedNeighborhood`.
impl<Item: MetricSpace<Impl> + Clone, Impl> BestCandidate<Item, Impl>
for RadiusBasedNeighborhood<Item, Impl>
{
    type Output = HashMap<u32,Item::Distance>;

    #[inline]
    fn consider(
        &mut self,
        _: &Item,
        distance: Item::Distance,
        candidate_index: usize,
        _: &Item::UserData,
    ) {
        // If the distance is lower than the bound we include the index
        // in the result.
        if distance <= self.max_distance {
            self.ids.insert(candidate_index as u32);
            self.distances.insert(candidate_index as u32, distance);
        }
    }

    #[inline]
    fn distance(&self) -> Item::Distance {
        self.max_distance
    }
    fn result(self, _: &Item::UserData) -> Self::Output {
        self.distances
    }
}

pub fn vantage_point_string_graph(input_list: &InputList, progress: bool) -> StringGraph {

    let bar2: Option<ProgressBar> = if progress {
        info!("Processing input list into nodes (progress bar may end early due to duplicate IDs)");
        Some(ProgressBar::new(input_list.strings.len() as u64))
    } else {
        None
    };

    let mut graph = GraphMap::<u32, u32, Undirected>::new();
    let mut string_to_node: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut node_to_string: HashMap<u32, Vec<u8>> = HashMap::new();
    let mut current_index = 0;

    input_list.strings.iter().for_each(|str| {
        if !string_to_node.contains_key(str) {
            graph.add_node(current_index);
            string_to_node.insert(str.clone(), current_index);
            node_to_string.insert(current_index, str.clone());
            current_index += 1;

        }
    });

    let vp = vpsearch::Tree::new(&input_list.strings.iter().map(|t| U8String{ u8str: t.clone() }).collect::<Vec<U8String>>());

    string_to_node.iter().enumerate().for_each(|(index,(x,n))| {
        //let expected = HashSet::new();
        let nearest = vp.find_nearest_custom(
            &U8String{u8str: x.clone()},
            &(),
            RadiusBasedNeighborhood::new(input_list.max_dist.clone() as u32),
        );
        nearest.iter().for_each(|(index,dist)| {
            graph.add_edge(*n, *index, *dist);
        });
        if index % 5000 == 0 && bar2.is_some() {
            bar2.as_ref().unwrap().inc(5000);
        }
    });


    StringGraph { graph, string_to_node, node_to_string, max_distance: input_list.max_dist as u64 }
}


pub fn max_set_distance(set: &Vec<Vec<u8>>) -> u64 {
    let mut max = 0;
    set.iter().for_each(|l1| {
        set.iter().for_each(|l2| {
            max = u64::max(string_distance_no_break(l1, l2, &l1.len()) as u64, max);
        });
    });
    max
}

/// A heuristic approach to splitting over-connected graphs: find the most balanced split that
/// lowers each subgraph below the over-connected limit (2 * max_distance). If successful it returns
/// those strings, else None
#[allow(dead_code)]
pub fn split_subgroup(string_graph: &mut StringGraph) -> Option<Vec<Vec<Vec<u8>>>> {
    if max_set_distance(&string_graph.string_to_node.keys().map(|k| k.clone()).collect::<Vec<Vec<u8>>>()) <= string_graph.max_distance * 2 {
        return None
    }
    let mut best_balance: f64 = 1.0;
    let node_count = string_graph.graph.node_count() as f64;
    let mut best_left = Vec::new();
    let mut best_right = Vec::new();
    if string_graph.graph.node_count() > 200 {
        info!("Processing connected components for {} nodes", string_graph.graph.node_count());
    }
    for x in string_graph.graph.all_edges() {
        let mut new_graph = string_graph.graph.clone();
        new_graph.remove_edge(x.0, x.1);

        let comps = get_connected_components(&StringGraph {
            graph: new_graph,
            string_to_node: string_graph.string_to_node.clone(),
            node_to_string: string_graph.node_to_string.clone(),
            max_distance: string_graph.max_distance,
        });
        if comps.len() == 2 {
            let balance = f64::abs((comps[0].len() as f64) - (comps[1].len() as f64)) / node_count;
            let left_set_max_dist = max_set_distance(&comps[0]);
            let right_set_max_dist = max_set_distance(&comps[1]);

            if balance < best_balance && left_set_max_dist < (string_graph.max_distance * 2) && right_set_max_dist < (string_graph.max_distance * 2) {
                best_balance = balance;
                best_left = comps[0].clone();
                best_right = comps[1].clone();
            }
        }
    }
    if best_left.len() > 0 {
        Some(vec![best_left, best_right])
    } else {
        None
    }
}

pub fn get_connected_components(string_graph: &StringGraph) -> Vec<Vec<Vec<u8>>> {
    let graph_nodes: Vec<Vec<u32>> = tarjan_scc(&string_graph.graph);
    graph_nodes.into_iter().map(|nodes| {
        nodes.into_iter().map(|node| string_graph.node_to_string.get(&node).unwrap().clone()).collect::<Vec<Vec<u8>>>()
    }).collect::<Vec<Vec<Vec<u8>>>>()
}

#[cfg(test)]
mod tests {

    use std::time::Instant;
    use rand::distributions::{Slice, Uniform};
    use triple_accel::levenshtein_exp;
    use crate::utils::base_utils::edit_distance;
    use super::*;
    use crate::rand::distributions::Distribution;
    use crate::rand::Rng;

    #[test]
    fn string_distance_test() {
        let str1 = vec![b'A', b'A', b'A', b'A'];
        let str2 = vec![b'A', b'A', b'A', b'T'];

        let str_dist = string_distance_no_break(&str1, &str2, &str1.len());
        assert_eq!(1, str_dist);

        let str1 = vec![b'A', b'A', b'A', b'A'];
        let str2 = vec![b'A', b'A', b'A', b'A'];

        let str_dist = string_distance_no_break(&str1, &str2, &str1.len());
        assert_eq!(0, str_dist);

        let str1 = vec![b'T', b'T', b'T', b'T'];
        let str2 = vec![b'A', b'A', b'A', b'A'];


        let str_dist = string_distance_no_break(&str1, &str2, &str1.len());
        assert_eq!(4, str_dist);
    }


    // Function to generate a random nucleotide sequence of a given length
    fn generate_random_nucleotide_sequence(length: usize) -> Vec<u8> {
        let mut rng = rand::thread_rng();
        let nucleotides = b"ACGT"; // Byte string of nucleotides
        let dist = Uniform::from(0..nucleotides.len());

        (0..length)
            .map(|_| nucleotides[dist.sample(&mut rng)])
            .collect()
    }

    fn mutate_sequence(sequence: &Vec<u8>, error_rate: f64) -> Vec<u8> {
        let mut rng = rand::thread_rng();
        let nucleotides = b"ACGT";
        let dist = Uniform::from(0..nucleotides.len());

        let mut resulting_sequence = Vec::new();
        for (position,nucleotide) in sequence.iter().enumerate() {
            // Decide if this nucleotide should be mutated, based on the error rate
            if rng.gen::<f64>() < error_rate {
                let mut new_nucleotide = *nucleotide;
                // Ensure the new nucleotide is different from the original
                while new_nucleotide == *nucleotide {
                    new_nucleotide = nucleotides[dist.sample(&mut rng)];
                }
                resulting_sequence.push(new_nucleotide);
            } else {
                resulting_sequence.push(*nucleotide);
            }
        }
        resulting_sequence
    }

    fn generate_random_kmers_with_error_rate() {
        let str1 = vec![b'A', b'A', b'A', b'A'];
        let str2 = vec![b'A', b'A', b'A', b'T'];

        let str_dist = string_distance_no_break(&str1, &str2, &str1.len());
        assert_eq!(1, str_dist);

        let str1 = vec![b'A', b'A', b'A', b'A'];
        let str2 = vec![b'A', b'A', b'A', b'A'];

        let str_dist = string_distance_no_break(&str1, &str2, &str1.len());
        assert_eq!(0, str_dist);

        let str1 = vec![b'T', b'T', b'T', b'T'];
        let str2 = vec![b'A', b'A', b'A', b'A'];


        let str_dist = string_distance_no_break(&str1, &str2, &str1.len());
        assert_eq!(4, str_dist);
    }

    fn create_random_string(length: &usize) -> Vec<u8> {
        let vowels = ['A', 'C', 'G', 'T'];
        let vowels_dist = Slice::new(&vowels).unwrap();
        let rng = rand::thread_rng();

    // build a string of 10 vowels
        let vowel_string: String = rng
            .sample_iter(&vowels_dist)
            .take(*length)
            .collect();
        vowel_string.into_bytes()
    }

    fn create_set_of_random_strings(length: &usize, pool_size: &usize) -> Vec<Vec<u8>> {
        (0..*pool_size).map(|_| create_random_string(length)).into_iter().collect()
    }

    #[test]
    fn test_graph_creation_comp() {
        let now = Instant::now();

        for i in (1000..10000).step_by(5000) {
            println!("I {} ",i);
            let vec_of_vecs = create_set_of_random_strings(&20, &i);
            println!("Vec size {} one {} two {}", vec_of_vecs.len(),
                     String::from_utf8(vec_of_vecs.get(0).unwrap().clone()).unwrap(),
                     String::from_utf8(vec_of_vecs.get(1).unwrap().clone()).unwrap());

            let collection = InputList { strings: vec_of_vecs, max_dist: 1 };
            let _graph = vantage_point_string_graph(&collection, false);
            println!("string_distance_break high/low {}", now.elapsed().as_millis());
        }
    }

    #[test]
    fn test_edit_distance() {
        let str1 = vec![b'A', b'C', b'G', b'T', b'A'];
        let str2 = vec![b'A', b'C', b'G', b'T', b'A'];
        assert_eq!(edit_distance(&str1, &str2), 0);
        let str2 = vec![b'T', b'C', b'G', b'T', b'A'];
        assert_eq!(edit_distance(&str1, &str2), 1);
        let str2 = vec![b'a', b'C', b'G', b'T', b'A'];
        assert_eq!(edit_distance(&str1, &str2), 0);
        let str2 = vec![b'R', b'C', b'G', b'T', b'A'];
        assert_eq!(edit_distance(&str1, &str2), 0);
    }



    fn aln_distance(st1: &Vec<FastaBase>, st2: &Vec<FastaBase>) -> f64 {
        levenshtein_exp(&&FastaBase::vec_u8(&st1), &FastaBase::vec_u8(&st1)) as f64
        
    }

    #[test]
    fn test_sift4_vs_string_dist() {
        let low_match = "AAAAAAAA";
        let low_match_fb = FastaBase::from_str(&low_match);
        let low_match_bytes = low_match.to_string().into_bytes();
        let high_match = "TTTTTTTT";
        let high_match_fb = FastaBase::from_str(&high_match);
        let high_match_bytes = high_match.to_string().into_bytes();
        let close_match = "AATAAAAA";
        let close_match_bytes = close_match.to_string().into_bytes();

        let iterations = 1000000;
        for _i in 0..4 {
            let now = Instant::now();
            for _x in 0..iterations {
                let _ = aln_distance(&low_match_fb, &high_match_fb);
            }
            println!("Aligned high/low {}", now.elapsed().as_millis());

        }

        let iterations = 1000000;
        for _i in 0..4 {
            let now = Instant::now();
            for _x in 0..iterations {
                let sift_dist = sift4::simple(&low_match, &high_match);
                assert_eq!(sift_dist, 8);
            }
            println!("Sift4 high/low {}", now.elapsed().as_millis());

        }
        for _i in 0..4 {

            let now = Instant::now();
            for _x in 0..iterations {
                let _ = sift4::simple(&low_match, &close_match);
            }
            println!("Sift4 close/low {}", now.elapsed().as_millis());

        }
        for _i in 0..4 {
            let now = Instant::now();
            for _x in 0..iterations {
                let _ = string_distance_no_break(&low_match_bytes, &high_match_bytes, &2);
            }
            println!("string_distance high/low {}", now.elapsed().as_millis());

        }
        for _i in 0..4 {

            let now = Instant::now();
            for _x in 0..iterations {
                let _ = string_distance_no_break(&low_match_bytes, &close_match_bytes, &2);
            }
            println!("string_distance close/low {}", now.elapsed().as_millis());

        }
        for _i in 0..4 {

            let now = Instant::now();
            for _x in 0..iterations {
                let _ = string_distance_break(&low_match_bytes, &high_match_bytes, &2);
            }
            println!("string_distance_break high/low {}", now.elapsed().as_millis());

        }
        for _i in 0..4 {

            let now = Instant::now();
            for _x in 0..iterations {
                let _ = string_distance_break(&low_match_bytes, &close_match_bytes, &2);
            }
            println!("string_distance_break close/low {}", now.elapsed().as_millis());
        }
    }

}