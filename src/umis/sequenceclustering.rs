use std::borrow::Borrow;
use std::cmp::Ordering::Less;
use std::collections::{BTreeSet, HashMap};
use std::convert::TryInto;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::str;

use flate2::bufread::GzDecoder;
use petgraph::dot::Dot;
use petgraph::prelude::*;
use rand::{Rng, seq};
use rand::distributions::Standard;
use rand::prelude::*;

use indicatif::ProgressBar;

use crate::fasta_comparisons::*;
use crate::umis::bronkerbosch::BronKerbosch;

// 0.8

pub struct InputList {
    pub strings: Vec<Vec<u8>>,
    pub max_dist: u32, // the minimum distance to consider for creating an edge
}

pub struct StringGraph {
    pub graph: GraphMap<u32, u32, Undirected>,
    pub string_to_node: HashMap<Vec<u8>, u32>,
    pub node_to_string: HashMap<u32, Vec<u8>>,
}


pub fn string_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> u32 {
    assert_eq!(str1.len(), str2.len());
    str1.iter().zip(str2.iter()).map(|(c1, c2)| if c1 == c2 { 0 } else { 1 }).sum()
}

pub struct KnownList {
    pub known_list: Vec<Vec<u8>>,
    pub known_list_map: HashMap<Vec<u8>, BestHits>,
    known_list_subset: HashMap<Vec<u8>, Vec<Vec<u8>>>,
    known_list_subset_key_size: usize,
}

fn validate_barcode(barcode: &Vec<u8>) -> bool {
    barcode.iter().filter(|b| !KNOWNBASES.contains_key(*b)).map(|n| n).collect::<Vec<_>>().len() == 0 as usize
}

pub fn get_reader(path: &str) -> Result<Box<dyn BufRead>, &'static str> {
    let file_type = path.split(".").collect::<Vec<&str>>().last().unwrap().clone();

    match file_type {
        "gz" => {
            let reader = Box::new(GzDecoder::new(BufReader::new(File::open(path).expect("Unable to open input known file"))));
            Ok(Box::new(BufReader::new(reader)))
        }
        _ => {
            Ok(Box::new(BufReader::new(File::open(path).expect("Unable to open known input file"))))
        }
    }
}

pub fn load_knownlist(knownlist_file: &String, starting_nmer_size: usize) -> KnownList {
    let mut existing_mapping = HashMap::new();
    let mut known_list_subset: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();
    let mut test_set = Vec::new();

    let mut raw_reader = get_reader(knownlist_file).unwrap();

    let mut cnt = 0;

    let mut btree = BTreeSet::new();
    println!("Adding known barcodes...");
    for line in raw_reader.lines() {
        let bytes = line.unwrap().as_bytes().to_vec();
        if validate_barcode(&bytes) {
            btree.insert(bytes);
        }
    }

    // now the barcodes are in order, use this to generate grouped keys
    let mut prefix: Option<Vec<u8>> = None;
    let mut container : Vec<Vec<u8>> = Vec::new();
    for bytes in &btree {
        if !prefix.is_some() {prefix = Some(bytes[0..starting_nmer_size].to_vec());}

        test_set.push(bytes.clone());

        existing_mapping.insert(bytes.clone(), BestHits { hits: vec![bytes.clone()], distance: 0 });

        let first_x = bytes.clone()[0..starting_nmer_size].to_vec();
        if edit_distance(&first_x,prefix.as_ref().unwrap()) > 0 {
            known_list_subset.insert(prefix.unwrap(), container.clone());
            container.clear();
            prefix = Some(first_x);
        }
        container.push(bytes.clone());
    }
    known_list_subset.insert(prefix.unwrap(), container.clone());
    KnownList { known_list: test_set, known_list_map: existing_mapping, known_list_subset, known_list_subset_key_size: starting_nmer_size }
}

pub struct BestHits {
    pub hits: Vec<Vec<u8>>,
    pub distance: usize,
}

pub fn edit_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> usize {
    assert_eq!(str1.len(), str2.len());

    let mut dist: usize = 0;
    for i in 0..str1.len() {
        if !((DEGENERATEBASES.get(&str1[i]).is_some() && DEGENERATEBASES.get(&str1[i]).unwrap().contains_key(&str2[i])) ||
            (DEGENERATEBASES.get(&str2[i]).is_some() && DEGENERATEBASES.get(&str2[i]).unwrap().contains_key(&str1[i]))) {
            dist += 1;
        }
    }
    dist
}

pub fn correct_to_known_list(barcode: &Vec<u8>, kl: &mut KnownList, max_distance: usize) -> BestHits {
    let mut hits = Vec::new();
    let mut distance = max_distance;
    if kl.known_list_map.contains_key(barcode) {
        hits.push(barcode.clone());
        distance = 0;
        BestHits { hits, distance }
    } else {
        let mut min_distance = max_distance;
        let barcode_subslice = &barcode[0..kl.known_list_subset_key_size].to_vec();

        for candidate_key in kl.known_list_subset.keys() {
            let key_dist = edit_distance(barcode_subslice, &candidate_key);

            if key_dist <= min_distance {
                let subset = kl.known_list_subset.get(candidate_key).unwrap();

                for full_candidate in subset {
                    let dist = edit_distance(&full_candidate, barcode);
                    if dist < min_distance {
                        hits.clear();
                        min_distance = dist;
                    }
                    if dist == min_distance {
                        hits.push(full_candidate.clone());
                    }
                }
            }
        }
        kl.known_list_map.insert(barcode.clone(), BestHits { hits: hits.clone(), distance });
        BestHits { hits, distance }
    }
}


pub fn input_list_to_graph(input_list: &InputList, compare: fn(&Vec<u8>, &Vec<u8>) -> u32, progress: bool) -> StringGraph {
    let mut graph = GraphMap::<u32, u32, Undirected>::new();
    let mut string_to_node: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut node_to_string: HashMap<u32, Vec<u8>> = HashMap::new();


    let mut current_index = 0;

    let bar2: Option<ProgressBar> = if (progress) {
        println!("Processing input list into nodes (progress bar may end early due to duplicate IDs)");
        Some(ProgressBar::new((input_list.strings.len() as u64)))
    } else {
        None
    };

    input_list.strings.iter().for_each(|str| {
        if !string_to_node.contains_key(str) {
            graph.add_node(current_index);
            string_to_node.insert(str.clone(), current_index);
            node_to_string.insert(current_index, str.clone());
            current_index += 1;
            if bar2.is_some() {
                bar2.as_ref().unwrap().inc(1);
            }
        }
    });

    let bar: Option<ProgressBar> = if progress {
        println!("processing barcode-barcode distances (this can take a long time)...");
        Some(ProgressBar::new((string_to_node.len() as u64 * string_to_node.len() as u64) / (2 as u64)))
    } else {
        None
    };

    string_to_node.clone().into_iter().for_each(|(key1, node1)| {
        string_to_node.clone().into_iter().for_each(|(key2, node2)| {
            if key1 != key2 && key1.cmp(&key2) == Less {
                let dist = compare(&key1, &key2);
                if dist <= input_list.max_dist {
                    graph.add_edge(node1, node2, dist);
                }
                if bar.is_some() {
                    bar.as_ref().unwrap().inc(1);
                }
            }
        });
    });

    StringGraph { graph, string_to_node, node_to_string }
}

pub fn process_cliques(string_graph: &StringGraph) -> BronKerbosch<u32, u32> {
    let mut bronker = BronKerbosch::new(string_graph.graph.clone());
    bronker.compute();
    println!("Discovered {} cliques in the data", bronker.max_cliques.len());
    bronker
}


pub fn generate_random_string(length: usize) -> Vec<u8> {
    let bases = vec![b'A', b'C', b'G', b'T'];
    let mut rng = rand::thread_rng();
    let mut results = Vec::new();
    for _i in 0..length { results.push(bases.choose(&mut rand::thread_rng()).unwrap().to_owned()) }
    results
}


pub fn create_one_off_errors(template: &Vec<u8>) -> Vec<Vec<u8>> {
    let mut ret = Vec::new();
    let bases = vec![b'A', b'C', b'G', b'T'];
    ret.push(template.clone());
    for i in 0..template.len() {
        for b in &bases {
            if template[i] != *b {
                let mut str = Vec::new();
                str.extend(template[0..i].to_vec());
                let b_vec = vec![b];
                str.extend(b_vec);
                str.extend(template[i + 1..template.len()].to_vec());
                ret.push(str);
            }
        }
    }
    ret
}

pub fn permute_random_string(length: usize, error_rate: f64, count: usize) -> Vec<Vec<u8>> {
    let mut ret = Vec::new();
    let base_str = generate_random_string(length);
    let bases = vec![b'A', b'C', b'G', b'T'];
    ret.push(base_str.clone());
    let mut rng = rand::thread_rng();
    for _i in 0..count {
        ret.push(base_str.iter().map(|c| {
            if rng.gen::<f64>() <= error_rate {
                bases.choose(&mut rand::thread_rng()).unwrap().to_owned()
            } else {
                c.clone()
            }
        }).into_iter().collect());
    }
    ret
}

pub fn generate_simulated_data(length: usize, base_count: usize, error_strings_per_string: usize, base_error_pct: f64) -> Vec<Vec<u8>> {
    let mut returned_vec = Vec::new();
    for _i in 0..base_count {
        returned_vec.append(&mut permute_random_string(length, base_error_pct, error_strings_per_string));
    }
    returned_vec
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn string_distance_test() {
        let str1 = vec![b'A', b'A', b'A', b'A'];
        let str2 = vec![b'A', b'A', b'A', b'T'];

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(1, str_dist);

        let str1 = vec![b'A', b'A', b'A', b'A'];
        let str2 = vec![b'A', b'A', b'A', b'A'];

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(0, str_dist);

        let str1 = vec![b'T', b'T', b'T', b'T'];
        let str2 = vec![b'A', b'A', b'A', b'A'];


        let str_dist = string_distance(&str1, &str2);
        assert_eq!(4, str_dist);
    }

    #[test]
    fn test_graph_creation() {
        let vec_of_vecs = vec![vec![b'A', b'A', b'A', b'A'], vec![b'C', b'C', b'C', b'C'], vec![b'G', b'A', b'A', b'A'], vec![b'C', b'A', b'A', b'A']];
        let collection = InputList { strings: vec_of_vecs, max_dist: 6 };
        let graph = input_list_to_graph(&collection, string_distance, false);
        //println!("{}", Dot::new(&graph.graph));
        assert_eq!(6, graph.graph.edge_count()); // paths are
    }


    // #[test] dont run for now
    fn test_four_set_count() {
        let test_set = generate_simulated_data(10, 10, 10, 0.1);
        let mut group1 = create_one_off_errors(&generate_random_string(10));
        println!("Group 1 {}", group1.len());
        group1.extend(create_one_off_errors(&generate_random_string(10)));
        group1.extend(create_one_off_errors(&generate_random_string(10)));
        group1.extend(create_one_off_errors(&generate_random_string(10)));
        println!("Group 1 {}", group1.len());
        let graph = input_list_to_graph(&InputList { strings: group1, max_dist: 4 }, string_distance, false);
        println!("Graphing!");
        process_cliques(&graph);
    }


    #[test]
    fn test_larger_clique_count() {
        let test_set = generate_simulated_data(10, 10, 10, 0.1);
        println!("TEST SIZE {}", test_set.len());
        let graph = input_list_to_graph(&InputList { strings: test_set, max_dist: 4 }, string_distance, false);
        process_cliques(&graph);
    }

    #[test]
    fn test_existing_single_cell() {
        let mut test_set = Vec::new();

        if let Ok(lines) = read_lines("test_data/just_sequences_500_16s.txt") {
            // Consumes the iterator, returns an (Optional) String
            for line in lines {
                if let Ok(ip) = line {
                    test_set.push(ip.as_bytes().to_vec());
                }
            }
        }
        println!("TEST SIZE {}", test_set.len());
        let graph = input_list_to_graph(&InputList { strings: test_set, max_dist: 1 }, string_distance, false);
        println!("Making clique");
        process_cliques(&graph);
    }

    #[test]
    fn test_validate_barcode() {
        let test_vec = vec![b'A', b'T'];
        let is_valid = validate_barcode(&test_vec);
        assert_eq!(true, is_valid);

        let test_vec = vec![b'a', b'T'];
        let is_valid = validate_barcode(&test_vec);
        assert_eq!(true, is_valid);

        let test_vec = vec![b'A', b'E'];
        let is_valid = validate_barcode(&test_vec);
        assert_eq!(false, is_valid);

        let test_vec = vec![b'a', b'e'];
        let is_valid = validate_barcode(&test_vec);
        assert_eq!(false, is_valid);
    }

    /*
    #[test]
    fn test_10x_whitelist() {
        let valid_list = load_knownlist(&"test_data/3M-february-2018.txt.gz".to_string());

    }*/

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

    fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
        where P: AsRef<Path>, {
        let file = File::open(filename)?;
        Ok(io::BufReader::new(file).lines())
    }
}