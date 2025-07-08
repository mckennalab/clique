use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use serde::{Serialize,Deserialize};
use std::collections::{BTreeMap};

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone, Copy)]
pub enum UMISortType {
    KnownTag,
    DegenerateTag,
}
#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum MergeStrategy {
    Align,
    Concatenate,
    ConcatenateBothForward,
}

impl SequenceLayout {
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
    pub fn from_yaml(yaml_file: &str) -> SequenceLayout {

        let mut file = File::open(yaml_file).unwrap_or_else(|_x | panic!("Unable to open YAML configuration file: {}",yaml_file));

        let mut yaml_contents = String::new();

        file.read_to_string(&mut yaml_contents)
            .unwrap_or_else(|_x | panic!("Unable to read contents of YAML configuration file: {}",&yaml_file));

        let mut deserialized_map: SequenceLayout = serde_yaml::from_str(&yaml_contents).expect("Unable to de-yaml your input file");

        for reference in deserialized_map.references.values_mut() {

            let mut ordering = reference.umi_configurations.values().map(|umi_config| {
                umi_config.order
            }).collect::<Vec<usize>>();

            ordering.sort_by_key(|a| *a);

            assert!(ordering.iter().enumerate().all(|(i, order)| {
                i == *order
            }), "The UMIConfigurations must have sequential order numbers, starting at 0");

            assert_eq!(reference.target_types.len(), reference.targets.len(), "Target sequences and target type lists must be the same length");

            reference.fill_and_validate_target_positions();
        }
        
        
        deserialized_map
    }

    ///
    /// Validate that the reference sequence contains all of the bases that we need to extract UMIs
    ///
    /// # Arguments
    ///    * ref_bases - the reference sequence, as a vector of bases
    ///
    pub fn validate_reference_sequence(ref_bases: &[u8], configurations: &BTreeMap<String,UMIConfiguration>) -> bool {
        let existing_bases = ref_bases.iter().map(|base| (char::from(*base), true)).collect::<BTreeMap<_, _>>();

        configurations.iter().map(|(_name,umi_config)| {
            existing_bases.contains_key(&umi_config.symbol)
        }).all(|x| x)
    }

}
#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum AlignedReadOrientation {
    Forward,
    Reverse,
    ReverseComplement,
    Unknown,
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum ReadPosition {
    Read1 {
        orientation: AlignedReadOrientation
    },
    Read2 {
        orientation: AlignedReadOrientation
    },
    Index1 {
        orientation: AlignedReadOrientation
    },
    Index2 {
        orientation: AlignedReadOrientation
    },
    Spacer {
        spacer_sequence: String
    },
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub enum UMIPadding {
    Left,
    Right,
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct UMIConfiguration {
    pub symbol: char,
    pub file: Option<String>,
    pub reverse_complement_sequences: Option<bool>,
    pub sort_type: UMISortType,
    pub length: usize,
    pub order: usize,
    pub pad: Option<UMIPadding>,
    pub max_distance: usize,
    pub maximum_subsequences: Option<usize>,
    pub max_gaps: Option<usize>,
    pub minimum_collapsing_difference: Option<f64>,
    pub levenshtein_distance: Option<bool>,
}


#[derive(Debug, PartialEq, Hash, Serialize, Deserialize, Clone, Eq)]
pub enum TargetType {
    Static,
    Cas9WT,
    Cas12AWT,
    Cas9ABE,
    Cas9CBE,
    Cas9ABECBE,
    Cas12ABE,
    Cas12CBE,
    Cas12ABECBE,
    Cas9Homing,
    Cas9ABEPalindrome,
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct ReferenceRecord {
    pub sequence: String,
    pub umi_configurations: BTreeMap<String,UMIConfiguration>,
    pub targets: Vec<String>,
    pub target_types: Vec<TargetType>,
    pub target_locations: Option<Vec<usize>>,
}

impl ReferenceRecord {

    pub fn fill_and_validate_target_positions(&mut self) {
        assert!(self.target_locations.is_none());
        let mut positions = Vec::new();

        self.targets.iter().for_each(|target|
            positions.push(match self.sequence.find(target.as_str()) {
                Some(x) => x,
                None => {panic!("Unable to find target {} in reference {}, please check your target sequences", target, self.sequence)}
            }));

        self.target_locations = Some(positions);
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct SequenceLayout {
    pub aligner: Option<String>,
    pub merge: Option<MergeStrategy>,
    pub reads: Vec<ReadPosition>,
    pub known_strand: bool,
    pub references: BTreeMap<String,ReferenceRecord>,
}

impl SequenceLayout {
    pub fn get_sorted_umi_configurations(&self, reference_name: &String) -> Vec<UMIConfiguration> {
        let reference = self.references.get(reference_name);
        match reference {
            None => {
                panic!("Unable to find reference {}",reference_name);
            }
            Some(ref_obj) => {
                let mut presorted = ref_obj.umi_configurations.values().cloned().collect::<Vec<UMIConfiguration>>();
                presorted.sort_by(|a,b| a.order.cmp(&b.order));
                presorted
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_yaml_readback() {
        let configuration =
            SequenceLayout::from_yaml(&String::from("test_data/test_layout.yaml"));
        assert!(configuration.references.contains_key("shorter_reference"));
        assert!(configuration.references.get("shorter_reference").unwrap().umi_configurations.contains_key("cell_id"));
        assert_eq!(configuration.references.get("shorter_reference").unwrap().umi_configurations.get("cell_id").unwrap().symbol,'*');
    }


    #[test]
    #[should_panic]
    fn test_basic_yaml_readback_invalid_ordering() {
        SequenceLayout::from_yaml(&String::from("test_data/test_layout_invalid.yaml"));
    }

    #[test]
    #[should_panic]
    fn test_basic_yaml_readback_invalid_ordering2() {
        SequenceLayout::from_yaml(&String::from("test_data/test_layout_invalid2.yaml"));
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
