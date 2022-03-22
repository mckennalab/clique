use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::str;

use rand::seq::SliceRandom;

lazy_static! {
    static ref OTHERBASES: HashMap<u8, Vec<u8>> = {
        let mut hashedvalues = HashMap::new();
        hashedvalues.insert(b'a', vec!(b'C',b'G',b'T'));
        hashedvalues.insert(b'A', vec!(b'C',b'G',b'T'));
        hashedvalues.insert(b'c', vec!(b'A',b'G',b'T'));
        hashedvalues.insert(b'C', vec!(b'A',b'G',b'T'));
        hashedvalues.insert(b'g', vec!(b'C',b'A',b'T'));
        hashedvalues.insert(b'G', vec!(b'C',b'A',b'T'));
        hashedvalues.insert(b't', vec!(b'C',b'G',b'A'));
        hashedvalues.insert(b'T', vec!(b'C',b'G',b'A'));
        hashedvalues
    };
}

pub struct KeyValueDualMapping {
    pub key_to_value: HashMap<String, String>,
    pub value_to_key: HashMap<String, String>,
}

// stolen from https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html
pub fn read_known_list_file(filename: &String) -> Option<KeyValueDualMapping> {
    // guess at the type -- single or double column

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(filename) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                let split = ip.split("\t");
                let vec: Vec<&str> = split.collect();
                match vec.len()  {
                    1 => return Some(read_value_file(filename)),
                    2 => return Some(read_key_value_file(filename)),
                    _ => panic!("Unknown number of columns: {}",vec.len()),
                }
            }
        }
    }
    None
}

// stolen from https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html
pub fn read_key_value_file(filename: &String) -> KeyValueDualMapping {
    let mut key_to_value = HashMap::new();
    let mut value_to_key = HashMap::new();
    print!("Read file {}", filename);

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(filename) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                let split = ip.split("\t");
                let vec: Vec<&str> = split.collect();
                if vec.len() != 2 {
                    println!("line: {:?} not in two parts", vec);
                }
                key_to_value.insert(vec[0].to_string(), vec[1].to_string());
                value_to_key.insert(vec[1].to_string(), vec[0].to_string());
            }
        }
    }
    KeyValueDualMapping { key_to_value: key_to_value, value_to_key: value_to_key }
}

pub fn read_value_file(filename: &String) -> KeyValueDualMapping {
    let mut key_to_value = HashMap::new();
    let mut value_to_key = HashMap::new();
    print!("Read file {}", filename);

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(filename) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                let split = ip.split("\t");
                let vec: Vec<&str> = split.collect();
                if vec.len() != 1 {
                    println!("line: {:?} not in two parts", vec);
                }
                key_to_value.insert(vec[0].to_string(), vec[0].to_string());
                value_to_key.insert(vec[0].to_string(), vec[0].to_string());
            }
        }
    }
    KeyValueDualMapping { key_to_value: key_to_value, value_to_key: value_to_key }
}



// choose an random base that's not the current base
#[inline]
pub fn base_utf8_to_random(base: u8) -> u8 {
    OTHERBASES[&base].choose(&mut rand::thread_rng()).unwrap().clone()
}

pub fn mutation_set(initial_value: String) -> HashSet<String> {
    let mut all_combinations = HashSet::new();
    for (index, base) in initial_value.as_bytes().iter().enumerate() {
        for alt_base in &OTHERBASES[&base] {
            let base_str = initial_value[..index].to_owned() + str::from_utf8(&[*alt_base]).unwrap();
            let full = base_str + &initial_value[(index+1)..].to_owned();
            assert_eq!(full.len(), initial_value.len(), "mutated string of the wrong size {}, should be the same as {}", full, initial_value);
            all_combinations.insert(full);
        }
    }
    all_combinations
}


// diversify the values of the key->value barcode mapping to include all the one-off error
// sequences from the known barcode list
pub fn create_known_errors(input_kv_mapping: KeyValueDualMapping) -> KeyValueDualMapping {
    // take each key/value mapping, mutate the value, and make a new mapping from this value to the key
    let mut new_value_to_key = HashMap::new();
    for (value, key) in input_kv_mapping.value_to_key {
        new_value_to_key.insert(value.clone(), key.clone());
        for one_off_error_str in mutation_set(value.clone()) {
            if new_value_to_key.contains_key(&one_off_error_str) {
                panic!("repeat {} -> {}, original {} -> {}",&one_off_error_str,key,&one_off_error_str,new_value_to_key[&one_off_error_str]);
            }
            new_value_to_key.insert(one_off_error_str, key.clone());
        }
    }
    KeyValueDualMapping { key_to_value: input_kv_mapping.key_to_value, value_to_key: new_value_to_key }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mutate_base() {
        let tenX = String::from("A");
        let mut_string = mutation_set(tenX);

        assert_eq!(3, mut_string.len());
        let mut all_combinations = HashSet::new();
        all_combinations.insert(String::from("C"));
        all_combinations.insert(String::from("G"));
        all_combinations.insert(String::from("T"));

        assert_eq!(mut_string,all_combinations);

    }

    #[test]
    fn test_mutate_two_bases() {
        let tenX = String::from("AA");
        let mut_string = mutation_set(tenX);

        assert_eq!(6, mut_string.len());
        let mut all_combinations = HashSet::new();
        all_combinations.insert(String::from("CA"));
        all_combinations.insert(String::from("GA"));
        all_combinations.insert(String::from("TA"));

        all_combinations.insert(String::from("AC"));
        all_combinations.insert(String::from("AG"));
        all_combinations.insert(String::from("AT"));

        assert_eq!(mut_string,all_combinations);

    }

    #[test]
    fn test_large_10x_values() {
        let tenX = String::from("data/100_barcode_test.txt");
        let key_value = read_key_value_file(&tenX);

        assert_eq!(100, key_value.key_to_value.len());
    }

    #[test]
    fn test_errors_10x_values() {
        let tenX = String::from("data/1_barcode_test.txt");
        let key_value = read_key_value_file(&tenX);

        let all_errors = create_known_errors(key_value);
        assert_eq!(1, all_errors.key_to_value.len());
        assert_eq!(49, all_errors.value_to_key.len());
    }

    #[test]
    fn test_errors_100_keys_GESTALT_values() {
        let tenX = String::from("data/gestalt_3_level_RT.txt");
        let key_value = read_value_file(&tenX);

        let all_errors = create_known_errors(key_value);
        assert_eq!(384, all_errors.key_to_value.len());
        assert_eq!(384 + (3 * 10 * 384), all_errors.value_to_key.len());
    }
}