use std::collections::HashMap;

use petgraph::prelude::*;
use std::cmp::Ordering::Less;
use crate::bronkerbosch::BronKerbosch;

pub struct InputList {
    pub strings: Vec<String>,
    pub max_dist: u32, // the minimum distance to consider for creating an edge
}

pub struct StringGraph {
    pub graph: GraphMap<u32, u32, Undirected>,
    pub string_to_node: HashMap<String, u32>,
}

pub struct Clique {
    pub consensus: String,
    pub mapping: HashMap<String, String>,
    pub max_dist: u32, // the minimum distance used for creating an edge
}

pub fn string_distance(str1: &String, str2: &String) -> u32 {
    assert_eq!(str1.len(), str2.len());
    str1.to_uppercase().chars().zip(str2.to_uppercase().chars()).map(|(c1, c2)| if c1 == c2 { 0 } else { 1 }).sum()
}

fn input_list_to_graph(input_list: &InputList, compare: fn(&String, &String) -> u32) -> StringGraph {
    let mut graph = GraphMap::<u32,u32,Undirected>::new();
    let mut string_to_node: HashMap<String, u32> = HashMap::new();

    let mut currentIndex = 0;
    input_list.strings.iter().for_each(|str| {
        println!("str {}",str);
        if !string_to_node.contains_key(str) {
            graph.add_node(currentIndex);
            string_to_node.insert(str.clone(), currentIndex);
            currentIndex += 1;
        }
    });

    string_to_node.clone().into_iter().for_each(|(key1, node1)| {
        string_to_node.clone().into_iter().for_each(|(key2, node2)| {
            if key1 != key2 && key1.cmp(&key2) == Less {
                let dist = compare(&key1, &key2);
                if dist <= input_list.max_dist {
                    println!("{} to {}, dist {}",key1, key2, dist);
                    graph.add_edge(node1, node2, dist);
                }
            }
        });
    });

    StringGraph { graph, string_to_node }
}

fn graph_to_clique(string_graph: &StringGraph) {
    let mut bronker = BronKerbosch::new(string_graph.graph.clone());
    bronker.compute();
    println!("CLIQWUES {}",bronker.max_cliques.len());
}
/*
#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;
    use rand::prelude::*;
    use rand::{seq, Rng}; // 0.8

    use std::borrow::Borrow;
    use rand::distributions::Standard;

    #[test]
    fn string_distance_test() {
        let str1 = String::from("ACGT");
        let str2 = String::from("AGGT");

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(1, str_dist);

        let str1 = String::from("AAAA");
        let str2 = String::from("CCCC");

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(4, str_dist);

        let str1 = String::from("AAAA");
        let str2 = String::from("aaaa");

        let str_dist = string_distance(&str1, &str2);
        assert_eq!(0, str_dist);
    }

    #[test]
    fn test_graph_creation() {
        let collection = InputList{strings: vec!["AAAA".to_string(),"CCCC".to_string(),"GGGG".to_string(),"TTTT".to_string()], max_dist: 6};
        let graph = input_list_to_graph(&collection,string_distance);
        //println!("{}", Dot::new(&graph.graph));
        assert_eq!(6, graph.graph.edge_count()); // paths are
    }

    fn generate_random_string(length: usize) -> String {
        let bases = vec!['A','C','G','T'];
        let mut rng = rand::thread_rng();
        let mut results = Vec::new();
        for _i in 0..length {results.push(bases.choose(&mut rand::thread_rng()).unwrap().to_owned())}
        results.into_iter().collect()
    }

    fn create_one_off_errors(length: usize) -> Vec<String> {
        let mut ret = Vec::new();
        let base_str = generate_random_string(length);
        let bases = vec!['A','C','G','T'];
        ret.push(base_str.clone());
        for i in 0..length {
            for b in bases {
                if base_str.chars()[i] != b {
                    //let front = base_str[0..i];
                    //let character =
                    //let back  = base_str[i..base_str.len()];

                    ret.push(base_str[0..i].append(vec![b]).append(base_str.chars()[i..base_str.len()]).into_iter().collect());
                }
            }
        }
        ret
    }

    fn permute_random_string(length: usize, error_rate: f64, count: usize) -> Vec<String> {
        let mut ret = Vec::new();
        let base_str = generate_random_string(length);
        let bases = vec!['A','C','G','T'];
        ret.push(base_str.clone());
        let mut rng = rand::thread_rng();
        for _i in 0..count {
            ret.push(base_str.chars().map(|c| {
                if rng.gen::<f64>() <= error_rate {
                    bases.choose(&mut rand::thread_rng()).unwrap().to_owned()
                } else {
                    c.clone()
                }
            }).collect::<String>())
        }
        ret
    }

    fn generate_simulated_data(length: usize, base_count: usize, error_strings_per_string: usize, base_error_pct: f64) -> Vec<String> {
        let mut returned_vec = Vec::new();
        for _i in 0..base_count {
            returned_vec.append(&mut permute_random_string(length,base_error_pct,error_strings_per_string));
        }
        returned_vec
    }

    #[test]
    fn test_clique_count() {
        let collection = InputList{strings: vec!["AAAA".to_string(),"ATAA".to_string(),"CTCC".to_string(),"CACC".to_string(),"CCCC".to_string(),"GGGG".to_string(),"TTTT".to_string()], max_dist: 2};
        let graph = input_list_to_graph(&collection,string_distance);
        println!("{}", Dot::new(&graph.graph));
        graph_to_clique(&graph);

    }

    #[test]
    fn test_larger_clique_count() {
        let test_set = generate_simulated_data(10, 10, 10, 0.1);
        println!("TEST SIZE {}",test_set.len());
        let graph = input_list_to_graph(&InputList{strings: test_set, max_dist: 4},string_distance);
        println!("{}", Dot::new(&graph.graph));
        graph_to_clique(&graph);

    }
}*/