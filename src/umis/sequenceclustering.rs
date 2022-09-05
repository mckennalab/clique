use std::collections::HashMap;

use petgraph::prelude::*;
use std::cmp::Ordering::Less;
use crate::umis::bronkerbosch::BronKerbosch;
use std::str;

pub struct InputList {
    pub strings: Vec<Vec<u8>>,
    pub max_dist: u32, // the minimum distance to consider for creating an edge
}

pub struct StringGraph {
    pub graph: GraphMap<u32, u32, Undirected>,
    pub string_to_node: HashMap<Vec<u8>, u32>,
    pub node_to_string: HashMap<u32,Vec<u8>>,
}

pub struct Clique {
    pub consensus: String,
    pub mapping: HashMap<Vec<u8>, Vec<u8>>,
    pub max_dist: u32, // the minimum distance used for creating an edge
}

pub fn string_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> u32 {
    assert_eq!(str1.len(), str2.len());
    str1.iter().zip(str2.iter()).map(|(c1, c2)| if c1 == c2 { 0 } else { 1 }).sum()
}

fn input_list_to_graph(input_list: &InputList, compare: fn(&Vec<u8>, &Vec<u8>) -> u32) -> StringGraph {
    let mut graph = GraphMap::<u32,u32,Undirected>::new();
    let mut string_to_node: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut node_to_string: HashMap<u32,Vec<u8>> = HashMap::new();
    let mut string_to_count: HashMap<Vec<u8>, u32> = HashMap::new();


    let mut currentIndex = 0;

    input_list.strings.iter().for_each(|str| {
        //println!("str {}",str::from_utf8(str).unwrap());
        if !string_to_node.contains_key(str) {
            graph.add_node(currentIndex);
            string_to_node.insert(str.clone(), currentIndex);
            node_to_string.insert(currentIndex,str.clone());
            currentIndex += 1;
        }
    });

    string_to_node.clone().into_iter().for_each(|(key1, node1)| {
        string_to_node.clone().into_iter().for_each(|(key2, node2)| {
            if key1 != key2 && key1.cmp(&key2) == Less {
                let dist = compare(&key1, &key2);
                if dist <= input_list.max_dist {
                    //println!("{} to {}, dist {}",str::from_utf8(&key1).unwrap(), str::from_utf8(&key2).unwrap(), dist);
                    graph.add_edge(node1, node2, dist);
                }
            }
        });
    });

    StringGraph { graph, string_to_node, node_to_string }
}

fn graph_to_clique(string_graph: &StringGraph) {
    let mut bronker = BronKerbosch::new(string_graph.graph.clone());
    bronker.compute();
    println!("CLIQWUES {}",bronker.max_cliques.len());
    let unknown = vec![b'A',b'T'];

    for clique in bronker.max_cliques {
        let tt: Vec<&str> = clique.iter().map(|k| str::from_utf8(string_graph.node_to_string.get(k).unwrap_or_else(||&unknown)).unwrap()).into_iter().collect();
        println!("Clique: {:?} {:?}",clique,tt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;
    use rand::prelude::*;
    use rand::{seq, Rng}; // 0.8

    use std::borrow::Borrow;
    use rand::distributions::Standard;
    use std::io;
    use std::fs::File;
    use std::io::BufRead;
    use std::path::Path;

    #[test]
    fn string_distance_test() {
        let str1 = vec![b'A',b'A',b'A',b'A'];
        let str2 = vec![b'A',b'A',b'A',b'T'];

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(1, str_dist);

        let str1 = vec![b'A',b'A',b'A',b'A'];
        let str2 = vec![b'A',b'A',b'A',b'A'];

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(0, str_dist);

        let str1 =vec![b'T',b'T',b'T',b'T'];
        let str2 = vec![b'A',b'A',b'A',b'A'];


        let str_dist = string_distance(&str1, &str2);
        assert_eq!(4, str_dist);
    }

    #[test]
    fn test_graph_creation() {
        let vec_of_vecs = vec![vec![b'A',b'A',b'A',b'A'],vec![b'C',b'C',b'C',b'C'],vec![b'G',b'A',b'A',b'A'],vec![b'C',b'A',b'A',b'A']];
        let collection = InputList{strings: vec_of_vecs, max_dist: 6};
        let graph = input_list_to_graph(&collection,string_distance);
        //println!("{}", Dot::new(&graph.graph));
        assert_eq!(6, graph.graph.edge_count()); // paths are
    }

    fn generate_random_string(length: usize) -> Vec<u8> {
        let bases = vec![b'A',b'C',b'G',b'T'];
        let mut rng = rand::thread_rng();
        let mut results = Vec::new();
        for _i in 0..length {results.push(bases.choose(&mut rand::thread_rng()).unwrap().to_owned())}
        results
    }

    fn create_one_off_errors(length: usize) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        let base_str = generate_random_string(length);
        let bases = vec![b'A',b'C',b'G',b'T'];
        ret.push(base_str.clone());
        for i in 0..length {
            for b in &bases {
                if base_str[i] != *b {
                    let mut str = Vec::new();
                    str.extend(base_str[0..i].to_vec());
                    let b_vec = vec![b];
                    str.extend(b_vec);
                    str.extend(base_str[i+1..base_str.len()].to_vec());
                    ret.push(str);
                }

            }
        }
        ret
    }

    fn permute_random_string(length: usize, error_rate: f64, count: usize) -> Vec<Vec<u8>> {
        let mut ret = Vec::new();
        let base_str = generate_random_string(length);
        let bases = vec![b'A',b'C',b'G',b'T'];
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

    fn generate_simulated_data(length: usize, base_count: usize, error_strings_per_string: usize, base_error_pct: f64) -> Vec<Vec<u8>> {
        let mut returned_vec = Vec::new();
        for _i in 0..base_count {
            returned_vec.append(&mut permute_random_string(length,base_error_pct,error_strings_per_string));
        }
        returned_vec
    }


    #[test]
    fn test_four_set_count() {
        let test_set = generate_simulated_data(10, 10, 10, 0.1);
        let mut group1 = create_one_off_errors(10);
        println!("Group 1 {}",group1.len());
        group1.extend(create_one_off_errors(10));
        group1.extend(create_one_off_errors(10));
        group1.extend(create_one_off_errors(10));
        println!("Group 1 {}",group1.len());
        let graph = input_list_to_graph(&InputList{strings: group1, max_dist: 4},string_distance);
        println!("Graphing!");
        graph_to_clique(&graph);

    }


    #[test]
    fn test_larger_clique_count() {
        let test_set = generate_simulated_data(10, 10, 10, 0.1);
        println!("TEST SIZE {}",test_set.len());
        let graph = input_list_to_graph(&InputList{strings: test_set, max_dist: 4},string_distance);
        graph_to_clique(&graph);

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
        println!("TEST SIZE {}",test_set.len());
        let graph = input_list_to_graph(&InputList{strings: test_set, max_dist: 1},string_distance);
        println!("Making clique");
        graph_to_clique(&graph);

    }

    fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
        where P: AsRef<Path>, {
        let file = File::open(filename)?;
        Ok(io::BufReader::new(file).lines())
    }
}