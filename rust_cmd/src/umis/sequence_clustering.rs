use std::cmp::Ordering::Less;
use std::collections::{HashMap};

use petgraph::algo::{tarjan_scc};
use petgraph::prelude::*;
use rand::{Rng};
use rand::prelude::*;
use utils::base_utils::simple_edit_distance;

use indicatif::ProgressBar;

use crate::umis::bronkerbosch::BronKerbosch;

pub struct InputList {
    pub strings: Vec<Vec<u8>>,
    pub max_dist: u64, // the minimum distance to consider for creating an edge
}

pub struct StringGraph {
    pub graph: GraphMap<u32, u32, Undirected>,
    pub string_to_node: HashMap<Vec<u8>, u32>,
    pub node_to_string: HashMap<u32, Vec<u8>>,
    pub max_distance: u64,
}


pub fn string_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> u64 {
    assert_eq!(str1.len(), str2.len());
    str1.iter().zip(str2.iter()).map(|(c1, c2)| if c1 == c2 { 0 } else { 1 }).sum()
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
            let key_dist = simple_edit_distance(barcode_subslice, &candidate_key);

            if key_dist <= min_distance {
                let subset = kl.known_list_subset.get(candidate_key).unwrap();

                for full_candidate in subset {
                    if full_candidate.len() == barcode.len() {
                        let dist = simple_edit_distance(&full_candidate, barcode);
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
        }
        kl.known_list_map.insert(barcode.clone(), BestHits { hits: hits.clone(), distance });
        BestHits { hits, distance }
    }
}
pub struct BestHits {
    pub hits: Vec<Vec<u8>>,
    pub distance: usize,
}

impl Clone for BestHits {
    fn clone(&self) -> BestHits {
        BestHits { hits: self.hits.clone(), distance: self.distance }
    }
}

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

pub fn input_list_to_graph(input_list: &InputList, compare: fn(&Vec<u8>, &Vec<u8>) -> u64, progress: bool) -> StringGraph {
    let mut graph = GraphMap::<u32, u32, Undirected>::new();
    let mut string_to_node: HashMap<Vec<u8>, u32> = HashMap::new();
    let mut node_to_string: HashMap<u32, Vec<u8>> = HashMap::new();

    let mut current_index = 0;

    let bar2: Option<ProgressBar> = if progress {
        trace!("Processing input list into nodes (progress bar may end early due to duplicate IDs)");
        Some(ProgressBar::new(input_list.strings.len() as u64))
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
        trace!("processing barcode-barcode distances (this can take a long time)...");
        Some(ProgressBar::new((string_to_node.len() as u64 * string_to_node.len() as u64) / (2 as u64)))
    } else {
        None
    };

    string_to_node.clone().into_iter().for_each(|(key1, node1)| {
        string_to_node.clone().into_iter().for_each(|(key2, node2)| {
            if key1 != key2 && key1.cmp(&key2) == Less {
                let dist = compare(&key1, &key2);
                if dist <= input_list.max_dist {
                    graph.add_edge(node1, node2, dist as u32);
                }
                if bar.is_some() {
                    bar.as_ref().unwrap().inc(1);
                }
            }
        });
    });

    StringGraph { graph, string_to_node, node_to_string, max_distance: input_list.max_dist as u64 }
}


pub fn process_cliques(string_graph: &StringGraph) -> BronKerbosch<u32, u32> {
    let mut bronker = BronKerbosch::new(string_graph.graph.clone());
    bronker.compute();
    trace!("Discovered {} cliques in the data", bronker.max_cliques.len());
    bronker
}


pub fn max_set_distance(set: &Vec<Vec<u8>>) -> u64 {
    let mut max = 0;
    set.iter().for_each(|l1| {
        set.iter().for_each(|l2| {
            max = u64::max(string_distance(l1, l2) as u64, max);
        });
    });
    max
}

pub fn average_set_distance(set: &Vec<Vec<u8>>) -> f64 {
    let mut dist = 0.0;
    set.iter().for_each(|l1| {
        set.iter().for_each(|l2| {
            dist += string_distance(l1, l2) as f64;
        });
    });
    dist / ((set.len() * set.len()) as f64)
}

pub fn interset_distance(left: &Vec<Vec<u8>>, right: &Vec<Vec<u8>>) -> f64 {
    let mut dist = 0.0;
    left.iter().for_each(|l1| {
        right.iter().for_each(|r1| {
            dist += string_distance(l1, r1) as f64;
        });
    });
    dist / ((left.len() * right.len()) as f64)
}

pub fn set_distances(left: &Vec<Vec<u8>>, right: &Vec<Vec<u8>>) -> (u64, u64, f64) {
    (max_set_distance(left), max_set_distance(right), interset_distance(left, right))
}

pub fn max_graph_set_differences(string_graph: &StringGraph) -> u64 {
    max_set_distance(&string_graph.node_to_string.values().map(|c| c.clone()).collect::<Vec<Vec<u8>>>())
}

/// A heuristic approach to splitting over-connected graphs: find the most balanced split that
/// lowers each subgraph below the over-connected limit (2 * max_distance). If successful it returns
/// those strings, else None
pub fn split_subgroup(string_graph: &mut StringGraph) -> Option<Vec<Vec<Vec<u8>>>> {
    if max_set_distance(&string_graph.string_to_node.keys().map(|k| k.clone()).collect::<Vec<Vec<u8>>>()) <= string_graph.max_distance * 2 {
        return None
    }
    let mut best_edge: Option<(u32, u32, &u32)> = None;
    let mut best_balance: f64 = 1.0;
    let mut best_left_right = (0, 0);
    let node_count = string_graph.graph.node_count() as f64;
    let mut best_left = Vec::new();
    let mut best_right = Vec::new();

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
                best_left_right = (comps[0].len(), comps[1].len());
                best_balance = balance;
                best_edge = Some(x);
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

pub fn generate_simulated_data(length: usize, groups: usize, error_strings_per_string: usize, base_error_pct: f64) -> Vec<Vec<u8>> {
    let mut returned_vec = Vec::new();
    for _i in 0..groups {
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

    #[test]
    fn test_cc() {
        let barcode_length = 8;

        let test_set = generate_simulated_data(barcode_length, 10, 10, 0.1);
        let mut group1 = create_one_off_errors(&generate_random_string(barcode_length));
        println!("Group 1 {}", group1.len());
        group1.extend(create_one_off_errors(&generate_random_string(barcode_length)));
        group1.extend(create_one_off_errors(&generate_random_string(barcode_length)));
        group1.extend(create_one_off_errors(&generate_random_string(barcode_length)));
        println!("Group 1 {}", group1.len());
        let graph = input_list_to_graph(&InputList { strings: group1, max_dist: 1 }, string_distance, false);
        println!("Graphing!");
        let cn = get_connected_components(&graph);
        println!("Graphing! {} ", cn.len());
        for group in cn {
            println!("group");
            for g in group {
                println!("Entry {}", String::from_utf8(g).unwrap());
            }
        }
    }

    #[test]
    fn test_cc_group_size() {
        let barcode_length = 8;
        let mut set1 = Vec::new();
        set1.push(vec![b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', ]);
        set1.push(vec![b'T', b'T', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', ]);
        set1.push(vec![b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'T', b'T', ]);
        set1.push(vec![b'A', b'A', b'A', b'A', b'A', b'A', b'T', b'T', b'T', b'T', ]);

        set1.push(vec![b'T', b'T', b'A', b'A', b'A', b'A', b'T', b'T', b'T', b'T', ]);
        set1.push(vec![b'T', b'T', b'A', b'A', b'T', b'T', b'T', b'T', b'T', b'T', ]);
        let mut graph = input_list_to_graph(&InputList { strings: set1, max_dist: 2 }, string_distance, false);
        assert!(!split_subgroup(&mut graph).is_some());

        let mut set2 = Vec::new();
        set2.push(vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C', ]);
        set2.push(vec![b'T', b'T', b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C', ]);
        set2.push(vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'T', b'T', ]);
        let mut graph = input_list_to_graph(&InputList { strings: set2, max_dist: 2 }, string_distance, false);
        assert!(!split_subgroup(&mut graph).is_some());

        let mut set3 = Vec::new();
        set3.push(vec![b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G', ]);
        set3.push(vec![b'T', b'T', b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G', ]);
        set3.push(vec![b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'T', b'T', ]);
        let mut graph = input_list_to_graph(&InputList { strings: set3, max_dist: 2 }, string_distance, false);
        assert!(!split_subgroup(&mut graph).is_some());
    }
}