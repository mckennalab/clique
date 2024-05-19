use std::cmp::min;
use rustc_hash::FxHashMap;

/// The nodes of a tree in our Starcode implementation
pub struct StarcodeNode {
    /// the nucleotide at this node in the tree
    base: u8,
    /// the mapping of nucleotides to daughter nodes
    daughters: FxHashMap<u8, StarcodeNode>,
    /// have we been traversed in a search?
    touched: bool,
    /// are we the terminal node representing a complete input sequence
    occupied: Option<String>,
}

pub struct SequenceCounts {
    sequence: String,
    count: u32,
}

const ALPHABET: [u8; 5] = [b'A', b'C', b'G', b'T', b'-'];
const SPECIAL_NULL: u8  = b'/';
impl StarcodeNode {
    pub fn new(base: u8) -> StarcodeNode {
        StarcodeNode{ base, daughters: Default::default(), touched: false, occupied: None }
    }

    pub fn root() -> StarcodeNode {
        StarcodeNode{ base: SPECIAL_NULL, daughters: Default::default(), touched: false, occupied: None }
    }
    pub fn create_daughters(&mut self) {
        for a in ALPHABET {
            assert!(!self.daughters.contains_key(&a));
            self.daughters.insert(a,StarcodeNode::new(a));
        }
    }

    pub fn insert_sequence(&mut self, original_seq: &[u8], offset: usize) {
        self.touched = true;
        if offset == original_seq.len() - 1 {
            assert!(self.daughters.is_empty()); // make sure we're at the bottom of the tree
            self.occupied = Some(String::from_utf8(original_seq.to_vec()).unwrap())
        } else {
            self.daughters.get_mut(&original_seq[offset]).unwrap().insert_sequence(original_seq,offset + 1);
        }
    }

    pub fn recursive_daughters(node: &mut StarcodeNode, current_depth: usize) {
        let mut daughters = FxHashMap::default();
        if current_depth > 0 {
            for a in ALPHABET {
                let mut daughter = StarcodeNode::new(a);
                StarcodeNode::recursive_daughters(&mut daughter, current_depth - 1);
                daughters.insert(a, daughter);
            }
            node.daughters = daughters;
        }
    }

    pub fn compute_score_vec(query_u8: &u8, reference_slice: &[u8]) -> Vec<usize> {
        reference_slice.iter().map(|x| if x == query_u8 {0} else {1}).collect()
    }

    pub fn create_tree_of_depth(depth: usize) -> StarcodeNode {
        assert!(depth <= 30,"We only allow trees of depth 30 or less");
        let mut root = StarcodeNode::root();
        StarcodeNode::recursive_daughters(&mut root, depth);
        root
    }

    pub fn search_at_depth(&mut self, depth: usize, query_position: &usize, query: &[u8], growing_ref: &[u8], max_dist: &usize) -> Vec<String> {
        if self.occupied.is_some() {
            let 
            vec!()
        }
        // we're at search depth
        if depth == 0 {
            self.daughters.iter_mut().map(|(base, daughter)| {
                let up = StarcodeNode::compute_score_vec(base, &query[0..query_position]);
                let left = StarcodeNode::compute_score_vec(&query[query_position], &growing_ref);
                let match_mismatch = if base == &query[query_position] {0} else {1};
                let score = match_mismatch + min(up[up.len()-1],left[left.len()-1]);

                if score <= max_dist {

                }

            })
        } else {
            self.daughters.iter_mut().map(|(_b, d)| d.search_at_depth(depth -1, remaining_query, max_dist)).collect()
        }
    }
}

struct BestMatches {
    sequences: Vec<String>,
    matches: FxHashMap<String,String>,
}

impl BestMatches {
    pub fn new(input_sequences: &Vec<String>) -> BestMatches {
        assert!(input_sequences.len() > 0);
        let mut sequences = input_sequences.clone();
        sequences.sort();
        let max_sequence = sequences.iter().map(|x| x.len()).max().unwrap();
        let mut tree = StarcodeNode::create_tree_of_depth(max_sequence);
        let mut last_sequence = &sequences[0].as_bytes();

        for sequence in sequences {
            let this_sequence = sequence.as_bytes();
            let shared_prefix = BestMatches::shared_prefix_length(last_sequence,this_sequence);

        }
        BestMatches{ sequences, matches: Default::default() }
    }

    pub fn shared_prefix_length(st1: &[u8], st2: &[u8]) -> usize {
        let mut prefix = 0;
        while prefix < st1.len() && st1[prefix] == st2[prefix] {
            prefix += 1;
        }
        prefix
    }
}


#[cfg(test)]
mod tests {
    use crate::umis::starcode_tree::StarcodeNode;

    #[test]
    pub fn create_empty_tree() {
        StarcodeNode::create_tree_of_depth(5);
    }


}

