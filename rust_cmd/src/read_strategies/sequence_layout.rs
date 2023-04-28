use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use serde::{Serialize,Deserialize};
use std::collections::BTreeMap;

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum UMISortType {
    /// This enum represents how we should extract molecular identifiers and known sequences from the
    /// stream of reads. It specifies the location (as aligned to the reference), the maximum distance
    /// you allow before not considering two sequences to be from the same source.
    KnownTag,
    DegenerateTag,
}

impl SequenceLayoutDesign {
    /// Load up a YAML document describing the layout of specific sequences within the reads. This configuration is specific
    /// to each sequencing platform and sequencing type (10X, sci, etc). The layout is described below.
    ///
    /// # Supported base tags
    ///
    /// *align*: (optional) * - contains one optional member, which can be _true_ or _false_
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
    pub fn from_yaml(yaml_file: String) -> Option<SequenceLayoutDesign> {

        let mut file = File::open(&yaml_file).expect(&format!("Unable to open YAML configuration file: {}",&yaml_file));

        let mut yaml_contents = String::new();

        file.read_to_string(&mut yaml_contents)
            .expect(&format!("Unable to read contents of YAML configuration file: {}",&yaml_file));

        let deserialized_map: SequenceLayoutDesign = serde_yaml::from_str(&yaml_contents).expect("Unable to de-yaml your input file");

        Some(deserialized_map)
    }
}
#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub enum ReadPosition {
    READ1,
    READ2,
    INDEX1,
    INDEX2
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct UMIConfiguration {
    pub symbol: char,
    pub file: Option<String>,
    pub sort_type: UMISortType,
    pub length: usize,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct SequenceLayoutDesign {
    pub reads: Vec<ReadPosition>,
    pub umi_configurations: BTreeMap<String,UMIConfiguration>,
}

// *********************************************       *********************************************
// ********************************************* TESTS *********************************************
// *********************************************       *********************************************

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_yaml_readback() {
        let configuration =
            SequenceLayoutDesign::from_yaml(String::from("test_data/test_layout.yaml")).unwrap();
        assert!(configuration.umi_configurations.contains_key("cell_id"));
        assert_eq!(configuration.umi_configurations.get("cell_id").unwrap().symbol,'*');
    }
}
