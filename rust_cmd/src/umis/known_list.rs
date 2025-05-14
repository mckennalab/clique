use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use log::info;

use serde::{Deserialize, Serialize};

use triple_accel::{levenshtein_exp};
use vpsearch::{Tree};
use vpsearch::{MetricSpace};
use utils::read_utils::reverse_complement;
use crate::read_strategies::sequence_layout::{UMIConfiguration};

use super::sequence_clustering::{RadiusBasedNeighborhood};

#[derive(Clone, Serialize, Deserialize, Hash, PartialEq, Eq)]
pub struct FastaString {
    pub fa: Vec<u8>,
    pub distance: u32,
    pub count: usize,
}

impl FastaString {
    pub fn new(fa: Vec<u8>) -> FastaString {
        FastaString { fa: fa.clone(), distance: u32::MAX, count: 0 }
    }
}

impl MetricSpace for FastaString {
    type UserData = ();
    //type UserData = HashMap<FastaString,bool>; // I dont want to deal with options this friday morning
    type Distance = u32;

    fn distance(&self, other: &Self, _consider: &()) -> Self::Distance {
        levenshtein_exp(&self.fa, &other.fa)// self.hamming_distance(other) //levenshtein_exp(&self.fa_u8, &other.fa_u8) //
    }
}

// #[derive(Serialize, Deserialize, Clone, )]
pub struct KnownList {
    vantage_tree: Tree<FastaString>,
    //, (), Owned<HashMap<FastaString,bool>>>,
    exact_matches: HashMap<FastaString, BestF32Hits>,
    input_list: Vec<FastaString>,
}

impl KnownList {
    pub fn new(umi_type: &UMIConfiguration, _starting_nmer_size: &usize) -> KnownList {
        let filename = umi_type.file.clone().unwrap();
        let filename = filename.as_str();

        info!("Setting up known list reader using a vantage tree for file {}; large files may take a long time", filename);

        let rev_comp = umi_type
            .reverse_complement_sequences
            .unwrap_or(false);

        let input_list = KnownList::create_input_set(filename, &rev_comp);

        //let vantage_tree = vpsearch::Tree::new_with_user_data_owned(&input_list, HashMap::new());
        let vantage_tree = vpsearch::Tree::new(&input_list);

        let exact_matches: HashMap<FastaString, BestF32Hits> =
            input_list.iter().map(|le|
                (le.clone(), BestF32Hits { hits: vec![le.fa.clone()], distance: 0 })).collect();

        KnownList {
            vantage_tree,
            exact_matches,
            input_list,
        }
    }

    pub fn create_input_set(filename: &str, reverse_comp: &bool) -> Vec<FastaString> {
        let raw_reader = BufReader::new(File::open(filename).expect(&format!("Unable to open input file {}",filename)));
        let mut input_set = Vec::new();
        for line in raw_reader.lines() {
            let mut bytes = line.unwrap().into_bytes();
            if *reverse_comp {
                bytes = reverse_complement(&bytes);
            }
            input_set.push(FastaString::new(bytes));
        }
        input_set
    }

    pub fn correct_to_known_list(
        &mut self,
        barcode: &Vec<u8>,
        max_distance: &u32,
    ) -> BestF32Hits {
        let string_rep = FastaString::new(barcode.clone());

        match self.exact_matches.get(&string_rep) {
            None => {
                let nearest = self.vantage_tree.find_nearest_custom(
                    &string_rep,
                    &(),
                    RadiusBasedNeighborhood::new(*max_distance),
                );

                let ret = BestF32Hits {
                    hits: nearest.iter().map(|(id, _dist)| self.input_list.get(*id as usize).unwrap().fa.clone()).collect(),
                    distance: nearest.iter().map(|(_id, dist)| *dist).next().unwrap_or(max_distance + 1),
                };

                self.exact_matches.insert(string_rep.clone(), ret.clone());
                ret
            }
            Some(x) => {
                x.clone()
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BestF32Hits {
    pub hits: Vec<Vec<u8>>,
    pub distance: u32,
}


#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };
    
    
    use crate::{
        read_strategies::sequence_layout::{UMIConfiguration, UMISortType},
        umis::known_list::{KnownList},
    };

    #[test]
    fn test_real_known_set() {
        let known_5p_list = UMIConfiguration {
            symbol: '0',
            file: Some("test_data/subset_barcode_list_500.txt".to_string()),
            reverse_complement_sequences: Some(false),
            sort_type: UMISortType::KnownTag,
            length: 16,
            order: 0,
            pad: None,
            max_distance: 0,
            maximum_subsequences: Some(25000),
            max_gaps: Some(1),
            minimum_collapsing_difference: None,
        };

        println!("Creating vantage point tree");
        let mut known_lookup = KnownList::new(&known_5p_list, &8);
        let file = File::open("test_data/subset_barcode_list_500.txt".to_string()).unwrap();
        let reader = BufReader::new(file);

        println!("Testing the top 100 events");
        for line in reader.lines().take(100) {
            //println!("{}", line?);
            assert_eq!(
                known_lookup.correct_to_known_list(
                    &line.unwrap().into_bytes(),
                    &1,
                )
                    .hits
                    .len(),
                1
            );
        }
        assert_eq!(
            known_lookup.correct_to_known_list(
                &"AAACCCAAGCAGATAA".as_bytes().to_vec(),
                &1,
            )
                .hits
                .len(),
            1
        );

        assert_eq!(
            known_lookup.correct_to_known_list(
                &"TAACCCAAGCAGATAT".as_bytes().to_vec(),
                &1,
            )
                .hits
                .len(),
            1
        );
    }
/*
    #[test]
    fn test_fastastring_edit_distance() {
        let str1 = FastaString::from_string(&"AAAAA".to_ascii_uppercase());
        let str2 = FastaString::from_string(&"AAAAA".to_ascii_uppercase());
        let consider = ();
        assert_eq!(str1.distance(&str2, &consider), 0);//, &HashMap::new()), 0);
        assert_eq!(str1.hamming_distance(&str2), 0);

        let str1 = FastaString::from_string(&"CTCCCCTTTCCC".to_ascii_uppercase());
        let str2 = FastaString::from_string(&"CCCCCCTTTCCC".to_ascii_uppercase());
        assert_eq!(str1.distance(&str2, &consider), 1);
        assert_eq!(str1.hamming_distance(&str2), 1);

        //CCCAAAGGGTTT CCCAAAGGGTTT
        let str1 = FastaString::from_string(&"TGTTTTTTTAAA".to_ascii_uppercase());
        let str2 = FastaString::from_string(&"TTTTTTTTTAAA".to_ascii_uppercase());
        assert_eq!(str1.distance(&str2, &consider), 1);
        assert_eq!(str2.distance(&str1, &consider), 1);
        assert_eq!(str1.hamming_distance(&str2), 1);

        let str1 = FastaString::from_string(&"AAAAA".to_ascii_uppercase());
        let str2 = FastaString::from_string(&"AAAAT".to_ascii_uppercase());
        assert_eq!(str1.distance(&str2, &consider), 1);
        assert_eq!(str1.hamming_distance(&str2), 1);

        let str1 = FastaString::from_string(&"AAAAA".to_ascii_uppercase());
        let str2 = FastaString::from_string(&"AAAA".to_ascii_uppercase());
        //assert_eq!(str1.distance(&str2, &consider), 1);
        //assert_eq!(str1.hamming_distance(&str2), 0);

        let str1 = FastaString::from_string(&"AAAAA".to_ascii_uppercase());
        let str2 = FastaString::from_string(&"TTTTT".to_ascii_uppercase());
        assert_eq!(str1.distance(&str2, &consider), 5);
        assert_eq!(str1.hamming_distance(&str2), 5);
    }*/
}
