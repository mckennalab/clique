//! Generate a clique read-structure YAML from an annotated GenBank file.
//!
//! # Convention
//!
//! A feature is "lineage-relevant" iff a configurable text tag (default
//! `lineage_target`) appears, case-insensitively, in its name — where the name
//! is the first present of the `/label`, `/gene`, `/note`, or `/product`
//! qualifiers, else the feature kind. Each selected feature becomes either a
//! CRISPR target or an extracted UMI/barcode:
//!
//! * **UMI/barcode** if it carries a `/clique_umi`, `/clique_symbol`, or
//!   `/clique_sort_type` qualifier, or its name contains `umi`, `barcode`, or
//!   `cell`. Its reference span is overwritten with a symbol character so the
//!   Rust extractor can pull it out. Parameters come from qualifiers with
//!   defaults: `/clique_symbol` (else auto-assigned from `0-9`),
//!   `/clique_sort_type` (`KnownTag`/`DegenerateTag`; defaults to `KnownTag`
//!   when `/clique_file` is set, else `DegenerateTag`), `/clique_max_distance`
//!   (default 1), `/clique_file`. `order` and `length` are derived (start
//!   order, feature span).
//! * **CRISPR target** otherwise. Its real bases become a `targets` entry and
//!   its `target_type` comes from `/clique_type` (else a type name found in the
//!   feature name, else `Cas9WT`). Targets may not overlap a UMI span.
//!   Prime-edit targets additionally require `/clique_edit_offset`,
//!   `/clique_ref`, and `/clique_alt`; `/clique_strand`,
//!   `/clique_call_flank`, `/clique_rtt`, and `/clique_scaffold` are optional.
//!
//! The emitted layout has a single reference (`--reference-name`, else the
//! GenBank LOCUS), `merge: ConcatenateBothForward`, `known_strand: true`, and a
//! single forward `Read1`; adjust those defaults by hand for paired/other runs.

use std::collections::{BTreeMap, BTreeSet};

use gb_io::seq::{Feature, Seq};

use crate::read_strategies::sequence_layout::{
    AlignedReadOrientation, MergeStrategy, PrimeEditSpec, ReadPosition, ReferenceRecord,
    SequenceLayout, TargetStrand, TargetType, UMIConfiguration, UMISortType,
};

/// Symbols auto-assigned to UMIs that do not declare a `/clique_symbol`.
const SYMBOL_POOL: &[char] = &[
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
];

pub struct GenbankToYamlOptions {
    /// Text that must appear in a feature's name for it to be included.
    pub tag: String,
    /// Reference name for the emitted layout; falls back to the GenBank LOCUS.
    pub reference_name: Option<String>,
}

impl Default for GenbankToYamlOptions {
    fn default() -> Self {
        GenbankToYamlOptions { tag: "lineage_target".to_string(), reference_name: None }
    }
}

/// The display name of a feature: first present of label/gene/note/product,
/// else the feature kind.
fn feature_name(feature: &Feature) -> String {
    for key in ["label", "gene", "note", "product", "standard_name"] {
        if let Some(v) = feature.qualifier_values(key).next() {
            if !v.is_empty() {
                return v.to_string();
            }
        }
    }
    feature.kind.to_string()
}

fn first_qualifier(feature: &Feature, key: &str) -> Option<String> {
    feature.qualifier_values(key).next().map(|s| s.to_string())
}

fn has_qualifier(feature: &Feature, key: &str) -> bool {
    feature.qualifier_values(key).next().is_some()
}

/// Parse a `TargetType` from a qualifier or name token (case/separator
/// insensitive).
pub fn parse_target_type(s: &str) -> Option<TargetType> {
    match s.to_lowercase().replace(['_', '-', ' '], "").as_str() {
        "cas9wt" => Some(TargetType::Cas9WT),
        "cas12awt" => Some(TargetType::Cas12AWT),
        "cas9abe" => Some(TargetType::Cas9ABE),
        "cas9cbe" => Some(TargetType::Cas9CBE),
        "cas9abecbe" => Some(TargetType::Cas9ABECBE),
        "cas12abe" => Some(TargetType::Cas12ABE),
        "cas12cbe" => Some(TargetType::Cas12CBE),
        "cas12abecbe" => Some(TargetType::Cas12ABECBE),
        "cas9homing" => Some(TargetType::Cas9Homing),
        "cas9abepalindrome" => Some(TargetType::Cas9ABEPalindrome),
        "primeedit" | "primeediting" => Some(TargetType::PrimeEdit),
        "static" => Some(TargetType::Static),
        _ => None,
    }
}

/// The known target-type tokens, longest first so `Cas9ABECBE` is matched
/// before `Cas9ABE` when scanning a feature name.
const TYPE_TOKENS: &[&str] = &[
    "Cas9ABEPalindrome",
    "PrimeEditing",
    "PrimeEdit",
    "Cas12ABECBE",
    "Cas9ABECBE",
    "Cas12ABE",
    "Cas12CBE",
    "Cas9Homing",
    "Cas12AWT",
    "Cas9ABE",
    "Cas9CBE",
    "Cas9WT",
    "Static",
];

fn target_type_from_name(name: &str) -> Option<TargetType> {
    let lower = name.to_lowercase();
    for token in TYPE_TOKENS {
        if lower.contains(&token.to_lowercase()) {
            return parse_target_type(token);
        }
    }
    None
}

fn parse_sort_type(s: &str) -> Option<UMISortType> {
    match s.to_lowercase().replace(['_', '-', ' '], "").as_str() {
        "knowntag" | "known" => Some(UMISortType::KnownTag),
        "degeneratetag" | "degenerate" => Some(UMISortType::DegenerateTag),
        _ => None,
    }
}

fn is_umi(feature: &Feature, name: &str) -> bool {
    if has_qualifier(feature, "clique_umi")
        || has_qualifier(feature, "clique_symbol")
        || has_qualifier(feature, "clique_sort_type")
    {
        return true;
    }
    let n = name.to_lowercase();
    ["umi", "barcode", "cell"].iter().any(|k| n.contains(k))
}

struct PendingUmi {
    start: usize,
    end: usize,
    name: String,
    symbol: Option<char>,
    sort_type: UMISortType,
    max_distance: usize,
    file: Option<String>,
}

struct PendingTarget {
    start: usize,
    end: usize,
    target_type: TargetType,
    prime_edit: Option<PrimeEditSpec>,
}

fn parse_prime_edit_spec(feature: &Feature, name: &str) -> Result<PrimeEditSpec, String> {
    let edit_offset = first_qualifier(feature, "clique_edit_offset")
        .ok_or_else(|| format!("prime-edit target '{}' requires /clique_edit_offset", name))?
        .parse::<usize>()
        .map_err(|_| format!("prime-edit target '{}' has an invalid /clique_edit_offset", name))?;
    let reference = first_qualifier(feature, "clique_ref")
        .ok_or_else(|| format!("prime-edit target '{}' requires /clique_ref", name))?;
    let alternate = first_qualifier(feature, "clique_alt")
        .ok_or_else(|| format!("prime-edit target '{}' requires /clique_alt", name))?;
    let strand = match first_qualifier(feature, "clique_strand")
        .unwrap_or_else(|| "Forward".to_string())
        .to_lowercase()
        .as_str()
    {
        "forward" | "+" => TargetStrand::Forward,
        "reverse" | "-" => TargetStrand::Reverse,
        value => {
            return Err(format!(
                "prime-edit target '{}' has unsupported /clique_strand '{}'",
                name, value
            ))
        }
    };
    let call_flank = match first_qualifier(feature, "clique_call_flank") {
        Some(value) => value.parse::<usize>().map_err(|_| {
            format!("prime-edit target '{}' has an invalid /clique_call_flank", name)
        })?,
        None => 10,
    };

    Ok(PrimeEditSpec {
        edit_offset,
        reference,
        alternate,
        strand,
        call_flank,
        rtt_sequence: first_qualifier(feature, "clique_rtt"),
        scaffold_sequence: first_qualifier(feature, "clique_scaffold"),
    })
}

/// Build a [`SequenceLayout`] from a GenBank record following the convention
/// documented on this module. Returns an error string describing the first
/// structural problem encountered.
pub fn genbank_to_layout(
    record: &Seq,
    opts: &GenbankToYamlOptions,
) -> Result<SequenceLayout, String> {
    let original: Vec<u8> = record.seq.to_ascii_uppercase();
    if original.is_empty() {
        return Err("GenBank record has no sequence.".to_string());
    }
    let tag_lower = opts.tag.to_lowercase();

    let mut umis: Vec<PendingUmi> = Vec::new();
    let mut targets: Vec<PendingTarget> = Vec::new();

    for feature in &record.features {
        let name = feature_name(feature);
        if !name.to_lowercase().contains(&tag_lower) {
            continue;
        }
        let (start, end) = feature
            .location
            .find_bounds()
            .map_err(|e| format!("feature '{}' has an unusable location: {:?}", name, e))?;
        // gb_io's `find_bounds` collapses a `complement(a..b)` location to the same numeric
        // bounds as the forward strand, so target sequences remain reference-oriented. Prime-edit
        // incorporation direction is supplied explicitly through /clique_strand.
        if start < 0 || end < 0 || (end as usize) > original.len() || start >= end {
            return Err(format!(
                "feature '{}' bounds {}..{} fall outside the sequence (len {})",
                name, start, end, original.len()
            ));
        }
        let (start, end) = (start as usize, end as usize);

        if is_umi(feature, &name) {
            let sort_type = first_qualifier(feature, "clique_sort_type")
                .and_then(|s| parse_sort_type(&s))
                .unwrap_or_else(|| {
                    if has_qualifier(feature, "clique_file") {
                        UMISortType::KnownTag
                    } else {
                        UMISortType::DegenerateTag
                    }
                });
            let symbol = first_qualifier(feature, "clique_symbol")
                .and_then(|s| s.chars().next());
            let max_distance = first_qualifier(feature, "clique_max_distance")
                .and_then(|s| s.parse::<usize>().ok())
                .unwrap_or(1);
            let file = first_qualifier(feature, "clique_file");
            umis.push(PendingUmi { start, end, name, symbol, sort_type, max_distance, file });
        } else {
            let target_type = first_qualifier(feature, "clique_type")
                .and_then(|s| parse_target_type(&s))
                .or_else(|| target_type_from_name(&name))
                .unwrap_or(TargetType::Cas9WT);
            let prime_edit = if target_type == TargetType::PrimeEdit {
                Some(parse_prime_edit_spec(feature, &name)?)
            } else {
                None
            };
            targets.push(PendingTarget { start, end, target_type, prime_edit });
        }
    }

    if umis.is_empty() && targets.is_empty() {
        return Err(format!(
            "No features contain the tag '{}'. Check the tag or your GenBank annotations.",
            opts.tag
        ));
    }

    // Deterministic ordering by start position.
    umis.sort_by_key(|u| u.start);
    targets.sort_by_key(|t| t.start);

    // Reject a target that overlaps a UMI span: symbol substitution would
    // corrupt its bases and the layout parser would fail to locate it.
    for t in &targets {
        for u in &umis {
            if t.start < u.end && u.start < t.end {
                return Err(format!(
                    "target at {}..{} overlaps UMI '{}' at {}..{}; targets and UMIs must be disjoint",
                    t.start, t.end, u.name, u.start, u.end
                ));
            }
        }
    }

    // Assign symbols, honoring explicit ones and drawing the rest from the pool
    // (excluding any symbol a feature claimed explicitly).
    let mut explicit = BTreeSet::new();
    for umi in &umis {
        if let Some(symbol) = umi.symbol {
            if !symbol.is_ascii_digit() {
                return Err(format!(
                    "UMI '{}' uses unsupported symbol '{}'; /clique_symbol must be an ASCII digit 0-9",
                    umi.name, symbol
                ));
            }
            if !explicit.insert(symbol) {
                return Err(format!(
                    "UMI '{}' reuses symbol '{}'; /clique_symbol values must be unique",
                    umi.name, symbol
                ));
            }
        }
    }
    let available: Vec<char> =
        SYMBOL_POOL.iter().copied().filter(|c| !explicit.contains(c)).collect();
    let mut pool_idx = 0usize;
    let mut umi_configs: BTreeMap<String, UMIConfiguration> = BTreeMap::new();
    let mut ref_bases = original.clone();

    for (order, u) in umis.iter().enumerate() {
        let symbol = match u.symbol {
            Some(c) => c,
            None => {
                if pool_idx >= available.len() {
                    return Err("ran out of UMI symbols to assign; at most 10 UMIs are supported per reference".to_string());
                }
                let c = available[pool_idx];
                pool_idx += 1;
                c
            }
        };
        let length = u.end - u.start;
        for i in u.start..u.end {
            ref_bases[i] = symbol as u8;
        }
        let config = UMIConfiguration {
            symbol,
            file: u.file.clone(),
            reverse_complement_sequences: None,
            sort_type: u.sort_type,
            length,
            order,
            pad: None,
            max_distance: u.max_distance,
            maximum_subsequences: None,
            max_gaps: None,
            minimum_collapsing_difference: None,
            levenshtein_distance: None,
        };
        // Use a unique config key derived from the feature name.
        let key = unique_key(&mut umi_configs, sanitize_key(&u.name));
        umi_configs.insert(key, config);
    }

    let sequence = String::from_utf8(ref_bases)
        .map_err(|_| "reference sequence is not valid UTF-8 after substitution".to_string())?;

    let target_strings: Vec<String> = targets
        .iter()
        .map(|t| String::from_utf8(original[t.start..t.end].to_vec()).unwrap())
        .collect();
    let target_types: Vec<TargetType> = targets.iter().map(|t| t.target_type.clone()).collect();
    let target_locations: Vec<usize> = targets.iter().map(|t| t.start).collect();
    let prime_edits = targets
        .iter()
        .enumerate()
        .filter_map(|(index, target)| target.prime_edit.clone().map(|spec| (index, spec)))
        .collect();

    let reference_name = opts
        .reference_name
        .clone()
        .or_else(|| record.name.clone())
        .unwrap_or_else(|| "reference".to_string());

    let mut record_out = ReferenceRecord {
        sequence,
        umi_configurations: umi_configs,
        targets: target_strings,
        target_types,
        target_locations: Some(target_locations),
        prime_edits,
    };
    record_out.fill_and_validate_target_positions();

    let mut references = BTreeMap::new();
    references.insert(reference_name, record_out);

    Ok(SequenceLayout {
        aligner: None,
        merge: Some(MergeStrategy::ConcatenateBothForward),
        reads: vec![ReadPosition::Read1 { orientation: AlignedReadOrientation::Forward }],
        known_strand: true,
        references,
    })
}

/// Turn a feature name into a YAML-friendly config key.
fn sanitize_key(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|c| if c.is_alphanumeric() { c.to_ascii_lowercase() } else { '_' })
        .collect();
    let trimmed = cleaned.trim_matches('_').to_string();
    if trimmed.is_empty() {
        "umi".to_string()
    } else {
        trimmed
    }
}

fn unique_key(existing: &BTreeMap<String, UMIConfiguration>, base: String) -> String {
    if !existing.contains_key(&base) {
        return base;
    }
    let mut i = 2;
    loop {
        let candidate = format!("{}_{}", base, i);
        if !existing.contains_key(&candidate) {
            return candidate;
        }
        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::reader::parse_slice;
    use std::io::Write;
    use tempfile::NamedTempFile;

    const GB: &str = "\
LOCUS       test_amplicon             60 bp    DNA     linear   SYN 07-JUL-2026
FEATURES             Location/Qualifiers
     misc_feature    1..16
                     /label=\"lineage_target cell barcode\"
                     /clique_sort_type=\"KnownTag\"
                     /clique_max_distance=\"2\"
     misc_feature    17..28
                     /label=\"lineage_target umi\"
                     /clique_sort_type=\"DegenerateTag\"
     misc_feature    30..52
                     /label=\"lineage_target site1\"
                     /clique_type=\"Cas9ABE\"
     misc_feature    5..10
                     /label=\"unrelated primer\"
ORIGIN
        1 aaaaaaaaaa aaaaaaggcc tgtcatctta gctaagatga caggtttttt tttttttccc
//
";

    const REPEATED_TARGET_GB: &str = "\
LOCUS       repeated_targets          16 bp    DNA     linear   SYN 07-JUL-2026
FEATURES             Location/Qualifiers
     misc_feature    1..4
                     /label=\"lineage_target first\"
                     /clique_type=\"Cas9WT\"
     misc_feature    9..12
                     /label=\"lineage_target second\"
                     /clique_type=\"Cas9WT\"
ORIGIN
        1 aaaaccccaa aacccc
//
";

    const PRIME_EDIT_GB: &str = "\
LOCUS       prime_edit_target         16 bp    DNA     linear   SYN 07-JUL-2026
FEATURES             Location/Qualifiers
     misc_feature    1..16
                     /label=\"lineage_target prime edit\"
                     /clique_type=\"PrimeEdit\"
                     /clique_edit_offset=\"6\"
                     /clique_ref=\"CC\"
                     /clique_alt=\"TT\"
                     /clique_strand=\"Reverse\"
                     /clique_call_flank=\"3\"
                     /clique_rtt=\"TTGG\"
                     /clique_scaffold=\"AACCGG\"
ORIGIN
        1 aaaaccccgg ggtttt
//
";

    fn parse_one() -> Seq {
        parse_slice(GB.as_bytes()).unwrap().into_iter().next().unwrap()
    }

    #[test]
    fn test_generates_umis_and_targets() {
        let layout =
            genbank_to_layout(&parse_one(), &GenbankToYamlOptions::default()).unwrap();
        assert_eq!(layout.references.len(), 1);
        let rec = layout.references.values().next().unwrap();

        // Two UMIs (barcode + umi), one target (Cas9ABE).
        assert_eq!(rec.umi_configurations.len(), 2);
        assert_eq!(rec.targets.len(), 1);
        assert_eq!(rec.target_types, vec![TargetType::Cas9ABE]);
    }

    #[test]
    fn test_repeated_targets_keep_genbank_feature_locations() {
        let record = parse_slice(REPEATED_TARGET_GB.as_bytes())
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let layout = genbank_to_layout(&record, &GenbankToYamlOptions::default()).unwrap();
        let reference = layout.references.values().next().unwrap();

        assert_eq!(reference.targets, vec!["AAAA", "AAAA"]);
        assert_eq!(reference.target_locations, Some(vec![0, 8]));
    }

    #[test]
    fn test_prime_edit_qualifiers_generate_explicit_specification() {
        let record = parse_slice(PRIME_EDIT_GB.as_bytes())
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let layout = genbank_to_layout(&record, &GenbankToYamlOptions::default()).unwrap();
        let reference = layout.references.values().next().unwrap();

        assert_eq!(reference.target_types, vec![TargetType::PrimeEdit]);
        assert_eq!(
            reference.prime_edits.get(&0),
            Some(&PrimeEditSpec {
                edit_offset: 6,
                reference: "CC".to_string(),
                alternate: "TT".to_string(),
                strand: TargetStrand::Reverse,
                call_flank: 3,
                rtt_sequence: Some("TTGG".to_string()),
                scaffold_sequence: Some("AACCGG".to_string()),
            })
        );
    }

    #[test]
    fn test_umi_spans_are_symbol_substituted() {
        let layout =
            genbank_to_layout(&parse_one(), &GenbankToYamlOptions::default()).unwrap();
        let rec = layout.references.values().next().unwrap();

        // First UMI symbol auto-assigned '0', 16bp; second '1', 12bp.
        let cell = rec.umi_configurations.values().find(|c| c.order == 0).unwrap();
        let umi = rec.umi_configurations.values().find(|c| c.order == 1).unwrap();
        assert_eq!(cell.length, 16);
        assert_eq!(umi.length, 12);
        // The reference must contain each UMI symbol run.
        assert!(rec.sequence.contains(&cell.symbol.to_string().repeat(16)));
        assert!(rec.sequence.contains(&umi.symbol.to_string().repeat(12)));
    }

    #[test]
    fn test_target_is_findable_in_reference() {
        let layout =
            genbank_to_layout(&parse_one(), &GenbankToYamlOptions::default()).unwrap();
        let rec = layout.references.values().next().unwrap();
        // Every target must be a substring of the (symbol-substituted) reference,
        // which is exactly what the layout parser later asserts.
        for t in &rec.targets {
            assert!(rec.sequence.contains(t), "target {} not in reference", t);
        }
    }

    #[test]
    fn test_orders_are_sequential_from_zero() {
        let layout =
            genbank_to_layout(&parse_one(), &GenbankToYamlOptions::default()).unwrap();
        let rec = layout.references.values().next().unwrap();
        let mut orders: Vec<usize> = rec.umi_configurations.values().map(|c| c.order).collect();
        orders.sort();
        assert_eq!(orders, vec![0, 1]);
    }

    #[test]
    fn test_untagged_features_are_ignored() {
        // The "unrelated primer" feature must not appear anywhere.
        let layout =
            genbank_to_layout(&parse_one(), &GenbankToYamlOptions::default()).unwrap();
        let rec = layout.references.values().next().unwrap();
        assert_eq!(rec.umi_configurations.len() + rec.targets.len(), 3);
    }

    #[test]
    fn test_custom_tag() {
        let opts = GenbankToYamlOptions { tag: "no_such_tag".to_string(), reference_name: None };
        let err = genbank_to_layout(&parse_one(), &opts).unwrap_err();
        assert!(err.contains("No features contain the tag"));
    }

    #[test]
    fn test_roundtrips_through_serde_yaml() {
        let layout =
            genbank_to_layout(&parse_one(), &GenbankToYamlOptions::default()).unwrap();
        let yaml = serde_yaml::to_string(&layout).unwrap();
        let mut yaml_file = NamedTempFile::new().unwrap();
        yaml_file.write_all(yaml.as_bytes()).unwrap();
        yaml_file.flush().unwrap();
        let back = SequenceLayout::from_yaml(yaml_file.path().to_str().unwrap());
        assert_eq!(back.references.len(), 1);
        assert!(yaml.contains("target_locations"));
        assert_eq!(
            back.references.values().next().unwrap().target_locations,
            layout.references.values().next().unwrap().target_locations
        );
    }

    #[test]
    fn test_rejects_unsupported_explicit_umi_symbol() {
        let invalid = GB.replace(
            "/clique_sort_type=\"KnownTag\"",
            "/clique_symbol=\"*\"",
        );
        let record = parse_slice(invalid.as_bytes())
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let error = genbank_to_layout(&record, &GenbankToYamlOptions::default()).unwrap_err();

        assert!(error.contains("must be an ASCII digit 0-9"));
    }

    #[test]
    fn test_rejects_duplicate_explicit_umi_symbols() {
        let invalid = GB
            .replace(
                "/clique_sort_type=\"KnownTag\"",
                "/clique_symbol=\"0\"",
            )
            .replace(
                "/clique_sort_type=\"DegenerateTag\"",
                "/clique_symbol=\"0\"",
            );
        let record = parse_slice(invalid.as_bytes())
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let error = genbank_to_layout(&record, &GenbankToYamlOptions::default()).unwrap_err();

        assert!(error.contains("values must be unique"));
    }
}
