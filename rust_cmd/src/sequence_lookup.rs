use std::fs::File;
use std::io::{BufRead, BufReader};
use symspell::{AsciiStringStrategy, SymSpell, Verbosity};
use crate::read_strategies::sequence_layout::UMIConfiguration;
use std::convert::TryFrom;

pub struct KnownLookup {
    corrector: SymSpell<AsciiStringStrategy>,
}

pub struct FastaStrategy {

}

impl KnownLookup {
    pub fn from(input_configuration: &UMIConfiguration) -> KnownLookup {
        let mut symspell: SymSpell<AsciiStringStrategy> = SymSpell::default();

        match &input_configuration.file {
            None => {panic!("unable to load known sequence list from a UMI configuration without a set filename")}
            Some(x) => {
                let file = File::open(x).expect("known UMI file not found");
                let sr = BufReader::new(file);
                for (_i, line) in sr.lines().enumerate() {
                    let mut line_str = line.unwrap();
                    line_str.push_str(" 1");
                    symspell.load_dictionary_line(&line_str, 0, 1, " ");
                }
            }
        }

        KnownLookup{ corrector: symspell}
    }

    pub fn correct(&self, sequence: &String, max_distance: &usize, if_multiple_take_first: bool) -> Option<String> {
        let result = self.corrector.lookup(sequence, Verbosity::Top, i64::try_from(*max_distance).unwrap());
        match result.len() {
            0 => None,
            1 => {Some(result[0].term.clone())}
            _ => {
                if if_multiple_take_first {
                    Some(result[0].term.clone())
                } else {
                    None
                }
            }
        }

    }
}

#[cfg(test)]
mod tests {
    use crate::read_strategies::sequence_layout::UMISortType::KnownTag;
    use super::*;
    use std::time::Instant;

    #[test]
    fn test_100k_by_100k_lookup() {
        let configuration = UMIConfiguration{
            symbol: '#',
            file: Some(String::from("test_data/100K-february-2018.txt")),
            reverse_complement_sequences: None,
            sort_type: KnownTag,
            length: 16,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: None,
            max_gaps: Some(0),
            minimum_collapsing_difference: None,
        };

        println!("loading file...");
        let now = Instant::now();
        let kf = KnownLookup::from(&configuration);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
        println!("Searching...");
        let now = Instant::now();
        for _x in 0..1000000 {
            let result = kf.corrector.lookup("AAACCCAAGAACCCGG", Verbosity::Top, 2);
            assert_eq!(result.len(),1);
        }
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
        // we used the top 100K entries in the 3M sequence file, so there's no sequences starting with Ts
        let result = kf.corrector.lookup("TTTCCCAAGAACCCGG", Verbosity::Top, 2);
        assert_eq!(result.len(),0);

    }
    /*
    #[test]
    fn test_3M_by_100K_lookup() {
        let configuration = UMIConfiguration{
            symbol: '#',
            file: Some(String::from("test_data/3M-february-2018.txt")),
            sort_type: KnownTag,
            length: 16,
        };

        println!("loading file...");
        let now = Instant::now();

        let kf = KnownLookup::from(&configuration);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);

        println!("Searching...");
        let now = Instant::now();

        for x in 0..100000 {
            let result = kf.corrector.lookup("AAACCCAAGAACCCGG", Verbosity::Top, 2);
            assert_eq!(result.len(),1);
        }
        let elapsed = now.elapsed();

        // we used the top 100K entries in the 3M sequence file, so there's no sequences starting with Ts
        let result = kf.corrector.lookup("TTTCCCAAGAACCCGG", Verbosity::Top, 2);
        assert_eq!(result.len(),1);

    }
*/

    #[test]
    fn test_simple_exact_correction() {

        let configuration = UMIConfiguration{
            symbol: '#',
            file: Some(String::from("test_data/just_sequences_500.txt")),
            reverse_complement_sequences: None,
            sort_type: KnownTag,
            length: 16,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: None,
            max_gaps: Some(0),
            minimum_collapsing_difference: None,
        };

        let kf = KnownLookup::from(&configuration);

        // we used the top 100K entries in the 3M sequence file, so there's no sequences starting with Ts
        let result = kf.correct(&String::from("ATATCCTAGACCCTGGGTGCTCCTTAG"), &2, false);
        assert_eq!(result.unwrap(),String::from("ATATCCTAGACCCTGGGTGCTCCTTAG"));

        let result = kf.correct(&String::from("AAAAACTAGACCCTGGGTGCTCCTTAG"), &2, false);
        assert_eq!(result.is_some(),false);

        for _x in 0..5000 {
            let _result = kf.correct(&String::from("AAAAACTAGACCCTGGGTGCTCCTTAG"), &2, false);
        }


    }
}
