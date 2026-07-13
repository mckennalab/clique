//! The read-structure YAML model. [`SequenceLayout`] describes each reference,
//! its UMI/barcode configurations, and its CRISPR targets; loading it validates
//! that symbols, targets, and UMI orderings are internally consistent.

use std::fmt::Debug;
use std::fs::File;
use std::io::Read;
use serde::{Serialize,Deserialize};
use std::collections::{BTreeMap, BTreeSet};

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

        for (reference_name, reference) in deserialized_map.references.iter_mut() {

            SequenceLayout::validate_umi_symbols(&reference.umi_configurations)
                .unwrap_or_else(|error| {
                    panic!(
                        "Invalid UMI configuration for reference '{}': {}",
                        reference_name, error
                    )
                });

            // Collect UMI `order` values. UMI names are BTreeMap keys and thus already unique --
            // duplicate YAML keys are dropped by serde during deserialization, so a name-collision
            // check here could never fire and is omitted.
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
        if SequenceLayout::validate_umi_symbols(configurations).is_err() {
            return false;
        }

        configurations.values().all(|umi_config| {
            ref_bases
                .iter()
                .filter(|base| char::from(**base) == umi_config.symbol)
                .count()
                == umi_config.length
        })
    }

    /// UMI symbols are embedded in both the reference and the second byte of
    /// `eX`/`oX` SAM tags. ASCII digits are the only characters that are both
    /// unambiguous reference markers and valid in that SAM tag position.
    pub fn validate_umi_symbols(
        configurations: &BTreeMap<String, UMIConfiguration>,
    ) -> Result<(), String> {
        let mut seen = BTreeSet::new();

        for (name, configuration) in configurations {
            if !configuration.symbol.is_ascii_digit() {
                return Err(format!(
                    "UMI '{}' uses unsupported symbol '{}'; symbols must be unique ASCII digits 0-9",
                    name, configuration.symbol
                ));
            }
            if !seen.insert(configuration.symbol) {
                return Err(format!(
                    "UMI '{}' reuses symbol '{}'; symbols must be unique within a reference",
                    name, configuration.symbol
                ));
            }
        }

        Ok(())
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reverse_complement_sequences: Option<bool>,
    pub sort_type: UMISortType,
    pub length: usize,
    pub order: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pad: Option<UMIPadding>,
    pub max_distance: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub maximum_subsequences: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_gaps: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub minimum_collapsing_difference: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub levenshtein_distance: Option<bool>,
}

impl UMIConfiguration {
    /// Known-tag correction uses Levenshtein distance unless explicitly disabled.
    pub fn uses_levenshtein_distance(&self) -> bool {
        self.levenshtein_distance.unwrap_or(true)
    }
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
    PrimeEdit,
}

#[derive(Debug, PartialEq, Hash, Serialize, Deserialize, Clone, Copy, Eq)]
pub enum TargetStrand {
    Forward,
    Reverse,
}

impl Default for TargetStrand {
    fn default() -> Self {
        TargetStrand::Forward
    }
}

/// Programmed prime-edit allele for one target. Coordinates and alleles are
/// always expressed in forward-reference orientation; `strand` determines the
/// direction used to recognize partial incorporation.
#[derive(Debug, PartialEq, Hash, Serialize, Deserialize, Clone, Eq)]
pub struct PrimeEditSpec {
    /// 0-based offset from the target start to the first replaced reference base.
    pub edit_offset: usize,
    /// Reference allele replaced by the edit. Empty for a programmed insertion.
    pub reference: String,
    /// Programmed allele. Empty for a programmed deletion.
    pub alternate: String,
    #[serde(default)]
    pub strand: TargetStrand,
    /// Reference bases retained on each side when classifying the haplotype.
    #[serde(default = "PrimeEditSpec::default_call_flank")]
    pub call_flank: usize,
    /// Optional RTT sequence in forward-reference orientation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rtt_sequence: Option<String>,
    /// Optional scaffold sequence used to recognize scaffold incorporation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scaffold_sequence: Option<String>,
}

impl PrimeEditSpec {
    fn default_call_flank() -> usize {
        10
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct ReferenceRecord {
    pub sequence: String,
    pub umi_configurations: BTreeMap<String,UMIConfiguration>,
    pub targets: Vec<String>,
    pub target_types: Vec<TargetType>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_locations: Option<Vec<usize>>,
    /// Prime-edit specifications keyed by target index.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub prime_edits: BTreeMap<usize, PrimeEditSpec>,
}

impl ReferenceRecord {

    pub fn fill_and_validate_target_positions(&mut self) {
        assert_eq!(
            self.target_types.len(),
            self.targets.len(),
            "Target sequences and target type lists must be the same length"
        );

        if let Some(positions) = self.target_locations.as_ref() {
            assert_eq!(
                positions.len(),
                self.targets.len(),
                "Target locations and target sequence lists must be the same length"
            );

            for (index, (target, start)) in self.targets.iter().zip(positions).enumerate() {
                assert!(!target.is_empty(), "Target {} must not be empty", index);
                let end = start.checked_add(target.len()).unwrap_or_else(|| {
                    panic!("Target {} location overflows the reference", index)
                });
                assert_eq!(
                    self.sequence.as_bytes().get(*start..end),
                    Some(target.as_bytes()),
                    "Target '{}' at location {} does not match reference sequence",
                    target,
                    start
                );
            }
            self.validate_prime_edits();
            return;
        }

        let mut positions = Vec::with_capacity(self.targets.len());
        let mut next_search_start = BTreeMap::new();

        for target in &self.targets {
            assert!(!target.is_empty(), "Target sequences must not be empty");
            let search_start = *next_search_start.get(target).unwrap_or(&0);
            let position = self
                .sequence
                .as_bytes()
                .get(search_start..)
                .and_then(|suffix| {
                    suffix
                        .windows(target.len())
                        .position(|window| window == target.as_bytes())
                })
                .map(|relative| search_start + relative)
                .unwrap_or_else(|| {
                    panic!(
                        "Unable to find occurrence of target {} at or after position {} in reference {}, please specify target_locations",
                        target,
                        search_start,
                        self.sequence
                    )
                });

            positions.push(position);
            next_search_start.insert(target.clone(), position + 1);
        }

        for (target, search_start) in next_search_start {
            let additional_position = self
                .sequence
                .as_bytes()
                .get(search_start..)
                .and_then(|suffix| {
                    suffix
                        .windows(target.len())
                        .position(|window| window == target.as_bytes())
                })
                .map(|relative| search_start + relative);

            if let Some(position) = additional_position {
                panic!(
                    "Target '{}' has an additional occurrence at position {}; please specify target_locations",
                    target,
                    position
                );
            }
        }

        self.target_locations = Some(positions);
        self.validate_prime_edits();
    }

    fn validate_prime_edits(&self) {
        let locations = self
            .target_locations
            .as_ref()
            .expect("Target locations must be resolved before validating prime edits");

        for target_index in self.prime_edits.keys() {
            assert!(
                *target_index < self.targets.len(),
                "Prime-edit specification references missing target index {}",
                target_index
            );
            assert_eq!(
                self.target_types[*target_index],
                TargetType::PrimeEdit,
                "Prime-edit specification at target index {} is attached to a non-prime target",
                target_index
            );
        }

        for (target_index, target_type) in self.target_types.iter().enumerate() {
            if target_type != &TargetType::PrimeEdit {
                continue;
            }
            let spec = self.prime_edits.get(&target_index).unwrap_or_else(|| {
                panic!(
                    "PrimeEdit target index {} must define a prime_edits entry",
                    target_index
                )
            });
            assert!(spec.call_flank > 0, "Prime-edit call_flank must be greater than zero");
            assert!(
                !(spec.reference.is_empty() && spec.alternate.is_empty()),
                "Prime edit at target index {} does not change the reference allele",
                target_index
            );
            validate_prime_sequence(&spec.reference, "reference allele");
            validate_prime_sequence(&spec.alternate, "alternate allele");
            if let Some(rtt) = spec.rtt_sequence.as_ref() {
                assert!(!rtt.is_empty(), "Prime-edit RTT sequence must not be empty");
                validate_prime_sequence(rtt, "RTT sequence");
            }
            if let Some(scaffold) = spec.scaffold_sequence.as_ref() {
                assert!(!scaffold.is_empty(), "Prime-edit scaffold sequence must not be empty");
                validate_prime_sequence(scaffold, "scaffold sequence");
            }

            let edit_start = locations[target_index]
                .checked_add(spec.edit_offset)
                .expect("Prime-edit start coordinate overflowed");
            let edit_end = edit_start
                .checked_add(spec.reference.len())
                .expect("Prime-edit end coordinate overflowed");
            let actual = self.sequence.as_bytes().get(edit_start..edit_end).unwrap_or_else(|| {
                panic!(
                    "Prime edit at target index {} falls outside the reference sequence",
                    target_index
                )
            });
            assert!(
                actual.eq_ignore_ascii_case(spec.reference.as_bytes()),
                "Prime-edit reference allele '{}' does not match reference sequence '{}' at position {}",
                spec.reference,
                String::from_utf8_lossy(actual),
                edit_start
            );
        }
    }
}

fn validate_prime_sequence(sequence: &str, description: &str) {
    assert!(
        sequence
            .bytes()
            .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N')),
        "Prime-edit {} contains a non-DNA base",
        description
    );
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct SequenceLayout {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligner: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
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
        assert_eq!(configuration.references.get("shorter_reference").unwrap().umi_configurations.get("cell_id").unwrap().symbol,'0');
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
    fn test_repeated_targets_infer_successive_occurrences() {
        let mut reference = ReferenceRecord {
            sequence: "AAAACCCCAAAA".to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec!["AAAA".to_string(), "AAAA".to_string()],
            target_types: vec![TargetType::Cas9WT, TargetType::Cas9WT],
            target_locations: None,
            prime_edits: BTreeMap::new(),
        };

        reference.fill_and_validate_target_positions();

        assert_eq!(reference.target_locations, Some(vec![0, 8]));
    }

    #[test]
    fn test_explicit_target_locations_are_preserved() {
        let mut reference = ReferenceRecord {
            sequence: "AAAACCCCAAAA".to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec!["AAAA".to_string(), "AAAA".to_string()],
            target_types: vec![TargetType::Cas9WT, TargetType::Cas9WT],
            target_locations: Some(vec![8, 0]),
            prime_edits: BTreeMap::new(),
        };

        reference.fill_and_validate_target_positions();

        assert_eq!(reference.target_locations, Some(vec![8, 0]));
    }

    #[test]
    #[should_panic(expected = "additional occurrence")]
    fn test_ambiguous_target_requires_explicit_location() {
        let mut reference = ReferenceRecord {
            sequence: "AAAACCCCAAAA".to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec!["AAAA".to_string()],
            target_types: vec![TargetType::Cas9WT],
            target_locations: None,
            prime_edits: BTreeMap::new(),
        };

        reference.fill_and_validate_target_positions();
    }

    #[test]
    #[should_panic(expected = "does not match reference sequence")]
    fn test_explicit_target_locations_must_match_reference() {
        let mut reference = ReferenceRecord {
            sequence: "AAAACCCCAAAA".to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec!["AAAA".to_string()],
            target_types: vec![TargetType::Cas9WT],
            target_locations: Some(vec![4]),
            prime_edits: BTreeMap::new(),
        };

        reference.fill_and_validate_target_positions();
    }

    fn configured_prime_edit(reference_allele: &str) -> ReferenceRecord {
        let mut prime_edits = BTreeMap::new();
        prime_edits.insert(0, PrimeEditSpec {
            edit_offset: 2,
            reference: reference_allele.to_string(),
            alternate: "TT".to_string(),
            strand: TargetStrand::Forward,
            call_flank: 2,
            rtt_sequence: Some("TTGG".to_string()),
            scaffold_sequence: Some("AACCGG".to_string()),
        });
        ReferenceRecord {
            sequence: "AACCGG".to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec!["AACCGG".to_string()],
            target_types: vec![TargetType::PrimeEdit],
            target_locations: Some(vec![0]),
            prime_edits,
        }
    }

    #[test]
    fn test_prime_edit_configuration_validates_and_roundtrips() {
        let mut reference = configured_prime_edit("CC");
        reference.fill_and_validate_target_positions();

        let yaml = serde_yaml::to_string(&reference).unwrap();
        let mut round_tripped: ReferenceRecord = serde_yaml::from_str(&yaml).unwrap();
        round_tripped.fill_and_validate_target_positions();
        assert_eq!(round_tripped, reference);
    }

    #[test]
    #[should_panic(expected = "does not match reference sequence")]
    fn test_prime_edit_rejects_wrong_reference_allele() {
        configured_prime_edit("GG").fill_and_validate_target_positions();
    }

    #[test]
    #[should_panic(expected = "must define a prime_edits entry")]
    fn test_prime_edit_requires_explicit_specification() {
        let mut reference = configured_prime_edit("CC");
        reference.prime_edits.clear();
        reference.fill_and_validate_target_positions();
    }


    #[test]
    fn test_validate_reference_sequence_all_present() {
        let mut configs = BTreeMap::new();
        configs.insert("umi1".to_string(), UMIConfiguration {
            symbol: '0',
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
        // Reference contains exactly the configured number of '0' markers.
        let ref_bases = b"ACGT0000000000ACGT";
        assert!(SequenceLayout::validate_reference_sequence(ref_bases, &configs));
    }

    #[test]
    fn test_validate_reference_sequence_missing_symbol() {
        let mut configs = BTreeMap::new();
        configs.insert("umi1".to_string(), UMIConfiguration {
            symbol: '1',
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
            symbol: '0',
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
            symbol: '1',
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
        // Has all ten '0' markers, but no '1' markers.
        let ref_bases = b"ACG0000000000TACGT";
        assert!(!SequenceLayout::validate_reference_sequence(ref_bases, &configs));

        // Has exactly ten '0' and five '1' markers.
        let ref_bases2 = b"ACG0000000000T11111ACGT";
        assert!(SequenceLayout::validate_reference_sequence(ref_bases2, &configs));
    }

    #[test]
    fn test_validate_reference_sequence_requires_exact_symbol_count() {
        let mut configs = BTreeMap::new();
        configs.insert("umi1".to_string(), UMIConfiguration {
            symbol: '0',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 3,
            order: 0,
            pad: None,
            max_distance: 1,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        });

        assert!(!SequenceLayout::validate_reference_sequence(b"AA00AA", &configs));
        assert!(SequenceLayout::validate_reference_sequence(b"AA000AA", &configs));
        assert!(!SequenceLayout::validate_reference_sequence(b"AA0000AA", &configs));
    }

    #[test]
    fn test_validate_umi_symbols_rejects_punctuation_and_duplicates() {
        let config = UMIConfiguration {
            symbol: '0',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 8,
            order: 0,
            pad: None,
            max_distance: 1,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        };
        let mut configs = BTreeMap::new();
        configs.insert("umi_1".to_string(), config.clone());
        assert!(SequenceLayout::validate_umi_symbols(&configs).is_ok());

        configs.get_mut("umi_1").unwrap().symbol = '*';
        assert!(SequenceLayout::validate_umi_symbols(&configs)
            .unwrap_err()
            .contains("unsupported symbol"));

        configs.get_mut("umi_1").unwrap().symbol = '0';
        let mut duplicate = config;
        duplicate.order = 1;
        configs.insert("umi_2".to_string(), duplicate);
        assert!(SequenceLayout::validate_umi_symbols(&configs)
            .unwrap_err()
            .contains("reuses symbol"));
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
            TargetType::Cas9Homing, TargetType::Cas9ABEPalindrome, TargetType::PrimeEdit,
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

    #[test]
    fn test_levenshtein_distance_defaults_to_true() {
        let mut config = UMIConfiguration {
            symbol: '0',
            file: Some("known_tags.txt".to_string()),
            reverse_complement_sequences: None,
            sort_type: UMISortType::KnownTag,
            length: 8,
            order: 0,
            pad: None,
            max_distance: 1,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        };

        assert!(config.uses_levenshtein_distance());
        config.levenshtein_distance = Some(true);
        assert!(config.uses_levenshtein_distance());
        config.levenshtein_distance = Some(false);
        assert!(!config.uses_levenshtein_distance());
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
