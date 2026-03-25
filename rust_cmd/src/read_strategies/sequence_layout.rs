use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use serde::{Serialize,Deserialize};
use std::collections::{BTreeMap, HashSet};

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

            let mut config_names : HashSet<String> = HashSet::default();
            
            let mut ordering = reference.umi_configurations.iter().map(|(name,umi_config)| {
                config_names.insert(name.clone());
                umi_config.order
            }).collect::<Vec<usize>>();

            assert_eq!(ordering.len(), config_names.len(), "Duplicate or mangled names in YAML configuration file; check umi_configurations field names and ordering");
            
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


    #[test]
    fn test_validate_reference_sequence_all_present() {
        let mut configs = BTreeMap::new();
        configs.insert("umi1".to_string(), UMIConfiguration {
            symbol: '*',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 10,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        });
        // Reference contains '*', so validation should pass
        let ref_bases = b"ACGT*ACGT";
        assert!(SequenceLayout::validate_reference_sequence(ref_bases, &configs));
    }

    #[test]
    fn test_validate_reference_sequence_missing_symbol() {
        let mut configs = BTreeMap::new();
        configs.insert("umi1".to_string(), UMIConfiguration {
            symbol: '#',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 10,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        });
        // Reference doesn't contain '#'
        let ref_bases = b"ACGTACGT";
        assert!(!SequenceLayout::validate_reference_sequence(ref_bases, &configs));
    }

    #[test]
    fn test_validate_reference_sequence_multiple_configs() {
        let mut configs = BTreeMap::new();
        configs.insert("umi1".to_string(), UMIConfiguration {
            symbol: '*',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 10,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        });
        configs.insert("umi2".to_string(), UMIConfiguration {
            symbol: '#',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::KnownTag,
            length: 5,
            order: 1,
            pad: None,
            max_distance: 1,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        });
        // Only has '*', not '#'
        let ref_bases = b"ACG*TACGT";
        assert!(!SequenceLayout::validate_reference_sequence(ref_bases, &configs));

        // Has both
        let ref_bases2 = b"ACG*T#ACGT";
        assert!(SequenceLayout::validate_reference_sequence(ref_bases2, &configs));
    }

    #[test]
    fn test_validate_reference_sequence_empty_configs() {
        let configs = BTreeMap::new();
        let ref_bases = b"ACGT";
        assert!(SequenceLayout::validate_reference_sequence(ref_bases, &configs));
    }

    #[test]
    fn test_umi_sort_type_serialization() {
        let known = UMISortType::KnownTag;
        let degen = UMISortType::DegenerateTag;
        let k_yaml = serde_yaml::to_string(&known).unwrap();
        let d_yaml = serde_yaml::to_string(&degen).unwrap();
        assert_ne!(k_yaml, d_yaml);
        let k_deser: UMISortType = serde_yaml::from_str(&k_yaml).unwrap();
        assert_eq!(k_deser, known);
    }

    #[test]
    fn test_merge_strategy_serialization() {
        let align = serde_yaml::to_string(&MergeStrategy::Align).unwrap();
        let concat = serde_yaml::to_string(&MergeStrategy::Concatenate).unwrap();
        let concat_fwd = serde_yaml::to_string(&MergeStrategy::ConcatenateBothForward).unwrap();
        // All should serialize differently
        assert_ne!(align, concat);
        assert_ne!(concat, concat_fwd);
    }

    #[test]
    fn test_aligned_read_orientation_variants() {
        assert_ne!(AlignedReadOrientation::Forward, AlignedReadOrientation::Reverse);
        assert_ne!(AlignedReadOrientation::Reverse, AlignedReadOrientation::ReverseComplement);
        assert_ne!(AlignedReadOrientation::ReverseComplement, AlignedReadOrientation::Unknown);
    }

    #[test]
    fn test_target_type_variants() {
        let types = vec![
            TargetType::Static, TargetType::Cas9WT, TargetType::Cas12AWT,
            TargetType::Cas9ABE, TargetType::Cas9CBE, TargetType::Cas9ABECBE,
            TargetType::Cas12ABE, TargetType::Cas12CBE, TargetType::Cas12ABECBE,
            TargetType::Cas9Homing, TargetType::Cas9ABEPalindrome,
        ];
        // All variants should be distinct
        for i in 0..types.len() {
            for j in (i + 1)..types.len() {
                assert_ne!(types[i], types[j]);
            }
        }
    }

    #[test]
    fn test_umi_padding_variants() {
        assert_ne!(UMIPadding::Left, UMIPadding::Right);
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
