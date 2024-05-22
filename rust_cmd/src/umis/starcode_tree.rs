use std::cmp::{max, min};
use std::collections::{HashMap, VecDeque};
use std::fmt;
use rustc_hash::FxHashMap;

use rand::prelude::*;

#[derive(Clone, Copy)]
pub struct NodeId {
    index: usize,
}

// Implement `Display` for `MinMax`.
impl fmt::Display for NodeId {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use `self.number` to refer to each positional data point.
        write!(f, "{}", self.index)
    }
}

type SequenceId = usize;


#[derive(Debug, Default)]
pub struct Associator {
    // node to node BY INDEX TO SAVE SPACE
    association_count: FxHashMap<SequenceId,usize>,
    associations: FxHashMap<SequenceId, Vec<SequenceId>>,
    next_id: SequenceId,
    underlying_sequences: FxHashMap<String, SequenceId>,
    underlying_index_to_string: FxHashMap<SequenceId, String>,
}

impl Associator {
    pub fn get_id(&mut self, id1: &String) -> SequenceId {
        match self.underlying_sequences.get(id1) {
            None => {
                self.underlying_sequences.insert(id1.clone(), self.next_id);
                self.underlying_index_to_string.insert(self.next_id, id1.clone());
                let ret = self.next_id.clone();
                self.next_id += 1;
                ret
            }
            Some(x) => { *x }
        }
    }
    /*
        pub fn cluster(&self) -> FxHashMap<usize, Vec<SequenceId>> {

        }
    */
    pub fn add_association(&mut self, id1: &String, id2: &String) -> (SequenceId, SequenceId) {
        let id1: SequenceId = self.get_id(id1);
        let id2: SequenceId = self.get_id(id2);
        self.add_association_to_ids(&id1, &id2)
    }
    pub fn add_association_to_id(&mut self, id1: &SequenceId, id2: &String) -> (SequenceId, SequenceId) {
        let id2: SequenceId = self.get_id(id2);
        self.add_association_to_ids(&id1, &id2)
    }
    pub fn add_association_to_ids(&mut self, id1: &SequenceId, id2: &SequenceId) -> (SequenceId, SequenceId) {
        {
            let id1_count = self.association_count.get(id1).unwrap_or(&0);
            self.association_count.insert(*id1, id1_count + 1);
        }
        {
            let id2_count = self.association_count.get(id2).unwrap_or(&0);
            self.association_count.insert(*id2, id2_count + 1);
        }

        let association = self.associations.entry(*id1).or_default();
        association.push(*id2);
        (*id1, *id2)
    }
}

struct TreeExplorePoint {
    node: NodeId,
    mismatch_count: usize,
    offset_into_seq: usize,
    depth: usize,
}

impl TreeExplorePoint {
    pub fn root() -> TreeExplorePoint {
        TreeExplorePoint { node: NodeId { index: 0 }, mismatch_count: 0, offset_into_seq: 0, depth: 0}
    }

    pub fn new(node: &NodeId, mismatch_count: &usize, offset_into_seq: &usize, depth: usize) -> TreeExplorePoint {
        TreeExplorePoint { node: *node, mismatch_count: *mismatch_count, offset_into_seq: *offset_into_seq, depth}
    }
}

pub struct Tree {
    tree_depth: usize,
    root: Node,
    max_mismatch: usize,
    associator: Associator,
    node_pool: NodePool,
}

pub struct NodePool {
    node_pool: Vec<Node>,
}

impl NodePool {
    pub fn new() -> NodePool {
        let root = NodePool::root();
        let mut pool = Vec::new();
        pool.push(root);
        NodePool { node_pool: pool }
    }

    pub fn create_remaining_node_string(&mut self, seq_database: &mut Associator, graft_node_index: &NodeId, graft_node_char_index: &usize, sequence: &[u8], offset: &usize) {
        let base_node = self.new_node(sequence[*offset]);
        let mut working_nodeid = base_node.clone();
        for i in (*offset + 1)..sequence.len() {
            let new_node = self.new_node(sequence[i]);
            match sequence[i] {
                b'A' => {
                    let mut working_node = &mut self.node_pool[working_nodeid.index];
                    working_node.acgt[0] = Some(new_node);
                }
                b'C' => {
                    let mut working_node = &mut self.node_pool[working_nodeid.index];
                    working_node.acgt[1] = Some(new_node);
                }
                b'G' => {
                    let mut working_node = &mut self.node_pool[working_nodeid.index];
                    working_node.acgt[2] = Some(new_node);
                }
                b'T' => {
                    let mut working_node = &mut self.node_pool[working_nodeid.index];
                    working_node.acgt[3] = Some(new_node);
                }
                _ => {}
            }
            working_nodeid = new_node.clone();
        }
        let str = String::from_utf8(sequence.to_vec()).unwrap();
        let seq_id = seq_database.get_id(&str);
        //println!("+++++++++++++++++++++         setting node {} to str id {}", working_nodeid,seq_id);
        self.node_pool.get_mut(working_nodeid.index).unwrap().sequence_offset = Some(seq_id);
        self.node_pool.get_mut(graft_node_index.index).unwrap().acgt[*graft_node_char_index] = Some(base_node);
    }

    pub fn new_node(&mut self, base: u8) -> NodeId {
        self.node_pool.push(Node { acgt: [None, None, None, None], empty_links: true, sequence_offset: None });
        NodeId { index: self.node_pool.len() - 1 }
    }

    fn root() -> Node {
        Node { acgt: [None, None, None, None], empty_links: true, sequence_offset: None }
    }
}

/// Prints the subtree starting from the node with the given `node_id`.
fn print_tree(arena: &NodePool, node_id: NodeId, depth: usize) {
    if let Some(node) = arena.node_pool.get(node_id.index) {
        let indent = "  ".repeat(depth);
        println!("{}Node Id: {}, Sequence Offset: {:?}", indent, node_id, node.sequence_offset);

        for (index, child_id) in node.acgt.iter().enumerate() {
            if let Some(child_id) = child_id {
                let nucleotide = match index {
                    0 => "A",
                    1 => "C",
                    2 => "G",
                    3 => "T",
                    _ => unreachable!(),
                };
                println!("{}{} -> Child Id: {}", indent, nucleotide, child_id);
                print_tree(arena, *child_id, depth + 1);
            }
        }
    }
}

impl Tree {
    pub fn new(depth: &usize, max_mismatch: &usize) -> Tree {
        Tree {
            tree_depth: *depth,
            root: NodePool::root(),
            max_mismatch: *max_mismatch,
            associator: Default::default(),
            node_pool: NodePool::new(),
        }
    }


    pub fn insert_sequence(&mut self,
                           sequence: &[u8]) {
        let mut nodes_to_explore: Vec<TreeExplorePoint> = vec![TreeExplorePoint::root()];
        let node_pool = &mut self.node_pool;
        let mut visited_nodes = 0;
        let mut asscoiated_nodes = 0;

        assert_eq!(sequence.len(), self.tree_depth);
        //println!("------------------------------------------------------------------- Called with {}------------------------------------------------------------------- ", String::from_utf8(sequence.to_vec()).unwrap());
        // we do two things:
        // 1) record any mismatch sequences hits in our associator database
        // 2) insert the input sequence into the tree at the corrent (no mismatch) location
        while !nodes_to_explore.is_empty() {
            visited_nodes += 1;
            let current_node = nodes_to_explore.pop().unwrap();
            //println!("Current node {} {} {} depth: {} -- {} - {} - {}", current_node.node.index, current_node.mismatch_count, current_node.offset_into_seq, current_node.depth,
            //         current_node.offset_into_seq == sequence.len(), current_node.mismatch_count <= self.max_mismatch, current_node.depth == self.tree_depth);

            if (current_node.offset_into_seq == sequence.len() && current_node.mismatch_count <= self.max_mismatch && current_node.depth == self.tree_depth) {
                //println!("adding");
                //if node_pool.node_pool.get(current_node.node.index).as_ref().unwrap().sequence_offset.is_none() {
                //println!(" -=-=-=-=- POINTS -=-=-=-=-=-=- {} -- {} {} -- {}",
                //         node_pool.node_pool.get(current_node.node.index).as_ref().unwrap().sequence_offset.is_none(),
                //         current_node.offset_into_seq,
                //         current_node.mismatch_count,
                //         current_node.node.index);

                //}
                let id2 = node_pool.node_pool.get(current_node.node.index).as_ref().unwrap().sequence_offset.unwrap();
                let idds = self.associator.add_association_to_id(&id2, &String::from_utf8(sequence.to_vec()).unwrap());
                asscoiated_nodes += 1;
            } else {
                let mut path_added = false;
                for i in 0..ALPHABETSIZE {
                    let x = node_pool.node_pool.get(current_node.node.index).unwrap().acgt[i];
                    let new_mismatch_count = current_node.mismatch_count + if current_node.offset_into_seq < sequence.len() && NODEALPHABET[i] == sequence[current_node.offset_into_seq] { 0 } else { 1 };

                    match (x,                                     // we already have the letter in our database
                           current_node.offset_into_seq < sequence.len(),     // we're still not at the end of our sequence
                           new_mismatch_count <= self.max_mismatch, // we have mismatches left to give
                           current_node.depth < self.tree_depth)              // we haven't exceeded the tree depth
                    {
                        (None, true, c, d) => {
                            //println!("A none false {} {}", c, d);
                            // we don't have the letter: this means we should make the path for the existing sequence, but that's the end of our adventure here -- there's nothing to match against
                            if NODEALPHABET[i] == sequence[current_node.offset_into_seq] && current_node.mismatch_count == 0 {
                                //println!("ADDDDDING");
                                node_pool.create_remaining_node_string(&mut self.associator, &current_node.node, &i, sequence, &current_node.offset_into_seq);
                            }
                            path_added = true;

                        },
                        (None, _, _, _) => {
                            //println!("B2 None, false, true, true");
                            // do nothing, we have no letters left to travese (tree paths to take, and
                        },
                        (Some(x), _, false, true) => {
                            //println!("B2a None, false, true, true");
                            // do nothing, we have no letters left to travese (tree paths to take, and
                        },
                        (Some(x), true, true, true) => {
                            //println!("C ");
                            // we have both sequence left and depth left - try match/mismatch, gap in reference, and gap in sequence
                            nodes_to_explore.push(TreeExplorePoint::new(&x, &new_mismatch_count, &(current_node.offset_into_seq + 1), current_node.depth + 1));
                            nodes_to_explore.push(TreeExplorePoint::new(&x, &(new_mismatch_count + 1), &current_node.offset_into_seq, current_node.depth + 1));
                            nodes_to_explore.push(TreeExplorePoint::new(&current_node.node, &(current_node.mismatch_count + 1), &(current_node.offset_into_seq + 1), current_node.depth));
                        },
                        (Some(x), true, true, false) => {
                            //println!("D ");
                            // we have sequence left but not depth left
                            nodes_to_explore.push(TreeExplorePoint::new(&current_node.node, &(current_node.mismatch_count + 1), &(current_node.offset_into_seq + 1), current_node.depth));
                        },
                        (Some(x), false, true, true) => {
                            //println!("E ");
                            // we have depth left but no sequence remaining
                            nodes_to_explore.push(TreeExplorePoint::new(&x, &(new_mismatch_count + 1), &current_node.offset_into_seq, current_node.depth + 1));
                        },
                        (a, b, c, d) => {
                            panic!("we're not handling {} {} {} {}", a.is_some(), b, c, d);
                        }
                    }
                }
            }
        }
        //println!("visited_nodes {} associated nodes {}",visited_nodes, asscoiated_nodes);
    }
}

const ALPHABETSIZE: usize = 4; // ensure that our node alphabet tracks our node internal storage
const NODEALPHABET: [u8; ALPHABETSIZE] = [b'A', b'C', b'G', b'T'];//, b'-'];
const SPECIAL_NULL: u8 = b'/';

/// The nodes of a tree in our Starcode implementation
pub struct Node {
    /// the mapping of nucleotides to daughter nodes
    acgt: [Option<NodeId>; ALPHABETSIZE],
    empty_links: bool,

    /// offset of this sequence in the underlying sorted array
    /// std::usize::MAX if not set
    sequence_offset: Option<SequenceId>,
}

impl Node {
    pub fn shared_prefix_length(str1: &[u8], str2: &[u8]) -> usize {
        let min = min(str1.len(), str2.len());
        for i in 0..min {
            if str1[i] != str2[i] {
                return i;
            }
        }
        min
    }
    /*
        pub fn create_tree(input_list: &mut Vec<String>) {
            assert!(input_list.len() > 1);

            /// pad sequences to the same length
            let max_length = input_list.iter().map(|x| x.len()).max();

            /// sort the list alphabetically
            input_list.sort();

            /// our base tree -- start with the first sequence
            let mut root = Node::root();
            root.insert_sequence(input_list[0].as_bytes(), 0, &0);

            let iter = input_list[..].windows(3);

            iter.for_each(|triplet| {
                let seed = Node::shared_prefix_length(triplet[1].as_bytes(),triplet[2].as_bytes());
                let start = Node::shared_prefix_length(triplet[0].as_bytes(),triplet[1].as_bytes());

                let mut hits : Vec<String> = Vec::new();

            });

            root.insert_sequence(input_list[0].as_bytes(),0,&0);


        }
    */
}

struct SearchThread<'a> {
    node: &'a Node,
    mismatches: usize,
    string_offset: usize,
}
/*
pub struct BestMatches {
    root: Node,
    matchings: HashMap<String, String>,
    max_distance: usize,
}

impl BestMatches {
    pub fn find_matches(&self, sequence: &[u8]) {
        let matches: Vec<String> = Vec::new();
        let mut search_points = VecDeque::from([SearchThread { node: &self.root, mismatches: 0, string_offset: 0 }]);

        while !search_points.is_empty() {
            let active_search_point = search_points.pop_front().unwrap();

            if active_search_point.node.daughters.is_empty() {}
            active_search_point.node.daughters.iter().for_each(|(base, node)| {
                if sequence[active_search_point.string_offset] == *base {}
            })
        }
    }
}*/


fn generate_combinations(length: usize) -> Vec<Vec<u8>> {
    // Define the nucleotides as bytes: A, C, G, T
    let nucleotides = vec![b'A', b'C', b'G', b'T'];

    // Recursive function to build combinations
    fn recurse(current: Vec<u8>, length: usize, nucleotides: &[u8], combinations: &mut Vec<Vec<u8>>) {
        if current.len() == length {
            combinations.push(current);
            return;
        }
        for &nucleotide in nucleotides {
            let mut new_combination = current.clone();
            new_combination.push(nucleotide);
            recurse(new_combination, length, nucleotides, combinations);
        }
    }

    let mut combinations = Vec::new();
    recurse(Vec::new(), length, &nucleotides, &mut combinations);
    combinations
}

#[cfg(test)]
mod tests {
    use rand::distributions::{Slice, Uniform};
    use rand::prelude::Distribution;
    use rand::Rng;

    use super::{Associator, generate_combinations, Node, NodeId, print_tree, SequenceId, Tree};

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
    fn mutate_sequence_set_times(sequence: &Vec<u8>, error_rate: f64) -> Vec<u8> {
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

    pub fn string_distance_no_break(str1: &Vec<u8>, str2: &Vec<u8>, _max_dist: &usize) -> usize {
        //assert_eq!(str1.len(), str2.len());
        str1.iter().zip(str2.iter()).map(|(c1, c2)| if c1 == c2 { 0 } else { 1 }).sum()
    }

    #[test]
    fn generate_random_kmers_with_error_rate() {
        let error_reads = 1000000;
        let sequence_size = 10;
        let str1 = vec![b'A'; sequence_size];
        let mut str1s : Vec<Vec<u8>> = (0..error_reads).map(|x| mutate_sequence(&str1, 0.1)).collect();
        let str2 = vec![b'C'; sequence_size];
        let mut str2s : Vec<Vec<u8>>  = (0..error_reads).map(|x| mutate_sequence(&str2, 0.1)).collect();
        let str3 = vec![b'G'; sequence_size];
        let mut str3s : Vec<Vec<u8>>  = (0..error_reads).map(|x| mutate_sequence(&str3, 0.1)).collect();
        let str4 = vec![b'T'; sequence_size];
        let mut str4s : Vec<Vec<u8>>  = (0..error_reads).map(|x| mutate_sequence(&str4, 0.1)).collect();
        str1s.append(&mut str2s);
        str1s.append(&mut str3s);
        str1s.append(&mut str4s);

        let mut tree = Tree::new(&sequence_size, &0);
        str1s.into_iter().enumerate().for_each(|combination| {
            tree.insert_sequence(combination.1.as_slice());

            if combination.0 % 10000 == 0 {
                println!("B count {}", combination.0);
            }
        });
        println!("tree size {}", tree.associator.associations.len());
        println!("tree size {}", tree.associator.association_count.len());
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
    pub fn test_basic_association() {
        let mut tree = Tree::new(&5, &1);
        tree.insert_sequence("ACTGA".as_bytes());
        print_tree(&tree.node_pool, NodeId { index: 0 }, 0);
        tree.insert_sequence("ACTGT".as_bytes());
        print_tree(&tree.node_pool, NodeId { index: 0 }, 0);

        assert_eq!(tree.associator.associations.len(), 1);
        /*
        -*/
    }

    #[test]
    pub fn test_insert_and_accumulate_sequences() {
        let sequence_size = 3;
        let all_combinations = generate_combinations(sequence_size);
        let mut tree = Tree::new(&sequence_size, &1);
        all_combinations.into_iter().enumerate().for_each(|combination| {
            tree.insert_sequence(combination.1.as_slice());

            if combination.0 % 10000 == 0 {
                println!("count {}", combination.0);
            }
        });
        println!("tree size {}", tree.associator.associations.len());
        println!("tree size {}", tree.associator.association_count.len());
    }


}

