use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use serde::{Serialize,Deserialize};
use std::collections::{BTreeMap};

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone, Copy)]
pub enum UMISortType {
    /// This enum represents how we should extract molecular identifiers and known sequences from the
    /// stream of reads. It specifies the location (as aligned to the reference), the maximum distance
    /// you allow before not considering two sequences to be from the same source.
    KnownTag,
    DegenerateTag,
}
#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum MergeStrategy {
    /// This enum represents how we should extract molecular identifiers and known sequences from the
    /// stream of reads. It specifies the location (as aligned to the reference), the maximum distance
    /// you allow before not considering two sequences to be from the same source.
    Align,
    Concatenate,
}


impl SequenceLayoutDesign {
    /// Load up a YAML document describing the layout of specific sequences within the reads. This configuration is specific
    /// to each sequencing platform and sequencing type (10X, sci, etc). The layout is described below.
    ///
    /// # Supported base tags
    ///
    /// *merge*: (optional) * - contains one optional member, which can be _align_ or _concatenate_
    /// *reads* - contains the read positions that are required for this configuration. The values, on individual lines, are _READ1, _READ2_, _INDEX1_, _INDEX2_
    /// *umi_configurations* - contains one section per UMI, each with:
    /// - *name* - the base of each UMI section
    ///   - *read* - which read we can extract this from
    ///   - *start* - the starting position, if align = true is set this this is in relation to the reference, otherwise the offset into the read
    ///   - *length* - how long this sequence is
    ///   - *file* - (optional) which file contains known sequences that we should match to. One sequence per line, no header
    ///
    /// an example of this format is the *test_layout.yaml* file in the test_data directory
    ///
    pub fn from_yaml(yaml_file: &String) -> Option<SequenceLayoutDesign> {

        let mut file = File::open(yaml_file).expect(&format!("Unable to open YAML configuration file: {}",yaml_file));

        let mut yaml_contents = String::new();

        file.read_to_string(&mut yaml_contents)
            .expect(&format!("Unable to read contents of YAML configuration file: {}",&yaml_file));

        let deserialized_map: SequenceLayoutDesign = serde_yaml::from_str(&yaml_contents).expect("Unable to de-yaml your input file");

        // validate that UMIConfigurations have orders and that they're sequential
        let mut ordering = deserialized_map.umi_configurations.iter().map(|(_name,umi_config)| {
            umi_config.order.clone()
        }).collect::<Vec<usize>>();
        ordering.sort_by(|a,b| a.cmp(b));

        assert!(ordering.iter().enumerate().all(|(i,order)| {
            i == *order
        }), "The UMIConfigurations must have sequential order numbers, starting at 0");

        Some(deserialized_map)
    }

    ///
    /// Validate that the reference sequence contains all of the bases that we need to extract UMIs
    ///
    /// # Arguments
    ///    * ref_bases - the reference sequence, as a vector of bases
    ///
    pub fn validate_reference_sequence(&self, ref_bases: &Vec<u8>) -> bool {
        let existing_bases = ref_bases.iter().map(|base| (char::from(base.clone()), true)).collect::<BTreeMap<_, _>>();

        self.umi_configurations.iter().map(|(_name,umi_config)| {
            existing_bases.contains_key(&umi_config.symbol)
        }).all(|x| x == true)
    }

}
#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum ReadPosition {
    READ1,
    READ2,
    INDEX1,
    INDEX2
}
impl ReadPosition {
    pub fn position(&self) -> usize {
        self.clone() as usize
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct UMIConfiguration {
    pub symbol: char,
    pub file: Option<String>,
    pub sort_type: UMISortType,
    pub length: usize,
    pub order: usize,
    pub max_distance: usize,
    pub maximum_subsequences: Option<usize>,
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct SequenceLayoutDesign {
    pub aligner: Option<String>,
    pub merge: Option<MergeStrategy>,
    pub reads: Vec<ReadPosition>,
    pub known_orientation: bool,
    pub read_separator: Option<String>,
    pub umi_configurations: BTreeMap<String,UMIConfiguration>,
}

impl SequenceLayoutDesign {
    pub fn get_sorted_umi_configurations(&self) -> Vec<UMIConfiguration> {
        let mut presorted = self.umi_configurations.iter().map(|(_name,umi_config)| {
            umi_config.clone()
        }).collect::<Vec<UMIConfiguration>>();
        presorted.sort_by(|a,b| a.order.cmp(&b.order));
        presorted
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_yaml_readback() {
        let configuration =
            SequenceLayoutDesign::from_yaml(&String::from("test_data/test_layout.yaml")).unwrap();
        assert!(configuration.umi_configurations.contains_key("cell_id"));
        assert_eq!(configuration.umi_configurations.get("cell_id").unwrap().symbol,'*');
    }


    #[test]
    #[should_panic]
    fn test_basic_yaml_readback_invalid_ordering() {
        let configuration =
            SequenceLayoutDesign::from_yaml(&String::from("test_data/test_layout_invalid.yaml")).unwrap();
        assert!(configuration.umi_configurations.contains_key("cell_id"));
        assert_eq!(configuration.umi_configurations.get("cell_id").unwrap().symbol,'*');
    }

    #[test]
    #[should_panic]
    fn test_basic_yaml_readback_invalid_ordering2() {
        let configuration =
            SequenceLayoutDesign::from_yaml(&String::from("test_data/test_layout_invalid2.yaml")).unwrap();
        assert!(configuration.umi_configurations.contains_key("cell_id"));
        assert_eq!(configuration.umi_configurations.get("cell_id").unwrap().symbol,'*');
    }


    /*
    TODO: figure out how to get SERDE to panic here or something else reasonable
    #[test]
    #[should_panic]
    fn test_parsing_wrong_data() {
        let configuration =
            SequenceLayoutDesign::from_yaml(String::from("test_data/test_layout-busted.yaml"));
    }
    */
}
