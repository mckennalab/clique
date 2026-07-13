//! Calling CRISPR edit "events" from an aligned read against its reference.
//!
//! This ports the intent of the (incomplete) Python `callers.py` into Rust.
//! Here we already have the gapped `reference_aligned` / `read_aligned`
//! strings and the per-reference target locations + types, so events are read
//! straight off the alignment columns rather than by re-parsing a CIGAR.
//!
//! An event is encoded in the McKenna-lab indel string format, all positions
//! 0-based in ungapped reference coordinates:
//!   * deletion:      `<len>D+<ref_pos>`
//!   * insertion:     `<len>I+<ref_pos>+<bases>`
//!   * substitution:  `<len>S+<ref_pos>+<bases>`   (base-editor scars)
//!   * `NONE` when a target site carries no called edit.
//!
//! Multiple events at one target are joined by `&`; the per-target calls for a
//! read are joined by `_` in target order (matching `Event.parse_event_string`
//! on the Python side, e.g. `10D+44_NONE_25D+76`).

use crate::read_strategies::sequence_layout::{
    PrimeEditSpec, ReferenceRecord, TargetStrand, TargetType,
};

const GAP: u8 = b'-';

/// The string emitted for a target site with no called edit.
pub const NONE_EVENT: &str = "NONE";

/// What kind of edit a target's chemistry produces, which decides which
/// alignment differences count as events inside the editing window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventClass {
    /// Nuclease double-strand break: insertions and deletions.
    DoubleStrandBreak,
    /// Adenine base editor: A->G (T->C on the opposite strand).
    AdenineBaseEdit,
    /// Cytosine base editor: C->T (G->A on the opposite strand).
    CytosineBaseEdit,
    /// Combined ABE+CBE: both substitution classes.
    DualBaseEdit,
    /// Programmed replacement evaluated against an explicit expected allele.
    PrimeEdit,
    /// No event calling (e.g. Static presence markers).
    None,
}

impl TargetType {
    /// Editing window as an inclusive `[start, end]` offset range within the
    /// target (0-based, relative to the target's start in the reference).
    ///
    /// Values follow the Python `callers.py` forward-strand windows. The
    /// Cas12 base-editor windows reuse the Cas9 base-editor window and should
    /// be recalibrated against real Cas12 base-editing data.
    pub fn editing_window(&self) -> (usize, usize) {
        match self {
            TargetType::Cas9WT | TargetType::Cas9Homing => (14, 19),
            TargetType::Cas12AWT => (14, 23),
            TargetType::Cas9ABE
            | TargetType::Cas9CBE
            | TargetType::Cas9ABECBE
            | TargetType::Cas9ABEPalindrome
            | TargetType::Cas12ABE
            | TargetType::Cas12CBE
            | TargetType::Cas12ABECBE => (2, 19),
            // Prime-edit windows come from PrimeEditSpec and are handled by
            // call_read_event_details rather than this legacy helper.
            TargetType::PrimeEdit => (0, 0),
            TargetType::Static => (0, 0),
        }
    }

    /// The class of edit this chemistry produces.
    pub fn event_class(&self) -> EventClass {
        match self {
            TargetType::Cas9WT | TargetType::Cas12AWT | TargetType::Cas9Homing => {
                EventClass::DoubleStrandBreak
            }
            TargetType::Cas9ABE | TargetType::Cas12ABE | TargetType::Cas9ABEPalindrome => {
                EventClass::AdenineBaseEdit
            }
            TargetType::Cas9CBE | TargetType::Cas12CBE => EventClass::CytosineBaseEdit,
            TargetType::Cas9ABECBE | TargetType::Cas12ABECBE => EventClass::DualBaseEdit,
            TargetType::PrimeEdit => EventClass::PrimeEdit,
            TargetType::Static => EventClass::None,
        }
    }
}

/// A single called edit relative to the reference. Positions are 0-based in
/// ungapped reference coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event {
    /// `length` reference bases deleted starting at `ref_pos`.
    Deletion { length: usize, ref_pos: usize },
    /// `bases` inserted immediately to the left of reference coordinate
    /// `ref_pos` (i.e. after `ref_pos` reference bases have been consumed).
    Insertion { ref_pos: usize, bases: Vec<u8> },
    /// A substituted base at `ref_pos`: reference had `ref_base`, read has
    /// `read_base`.
    Substitution { ref_pos: usize, ref_base: u8, read_base: u8 },
}

impl Event {
    /// Inclusive reference-coordinate span this event touches, for window
    /// overlap tests. An insertion is a zero-width point at `ref_pos`.
    fn ref_span(&self) -> (usize, usize) {
        match self {
            Event::Deletion { length, ref_pos } => {
                (*ref_pos, ref_pos + length.saturating_sub(1))
            }
            Event::Insertion { ref_pos, .. } => (*ref_pos, *ref_pos),
            Event::Substitution { ref_pos, .. } => (*ref_pos, *ref_pos),
        }
    }

    /// Render into the McKenna event string.
    pub fn encode(&self) -> String {
        match self {
            Event::Deletion { length, ref_pos } => format!("{}D+{}", length, ref_pos),
            Event::Insertion { ref_pos, bases } => {
                format!("{}I+{}+{}", bases.len(), ref_pos, String::from_utf8_lossy(bases))
            }
            Event::Substitution { ref_pos, read_base, .. } => {
                format!("1S+{}+{}", ref_pos, *read_base as char)
            }
        }
    }
}

/// Extract every indel and single-base substitution from a gapped alignment
/// pair. The two slices must be equal length. Substitutions are returned
/// unconditionally (including random mismatches); callers decide which to keep
/// based on the target chemistry and window.
pub fn extract_all_events(reference_aligned: &[u8], read_aligned: &[u8]) -> Vec<Event> {
    assert_eq!(
        reference_aligned.len(),
        read_aligned.len(),
        "aligned reference and read must be the same length"
    );

    let mut events = Vec::new();
    let mut ref_pos = 0usize; // next ungapped reference coordinate
    let mut i = 0usize;
    let n = reference_aligned.len();

    while i < n {
        let r = reference_aligned[i];
        let q = read_aligned[i];

        match (r == GAP, q == GAP) {
            // deletion run: reference present, read gap
            (false, true) => {
                let start = ref_pos;
                let mut length = 0;
                while i < n && read_aligned[i] == GAP && reference_aligned[i] != GAP {
                    length += 1;
                    ref_pos += 1;
                    i += 1;
                }
                events.push(Event::Deletion { length, ref_pos: start });
            }
            // insertion run: reference gap, read present
            (true, false) => {
                let start = ref_pos;
                let mut bases = Vec::new();
                while i < n && reference_aligned[i] == GAP && read_aligned[i] != GAP {
                    bases.push(read_aligned[i]);
                    i += 1;
                }
                events.push(Event::Insertion { ref_pos: start, bases });
            }
            // gap-gap column: nothing consumed on either side, skip
            (true, true) => {
                i += 1;
            }
            // aligned column: substitution if the bases differ
            (false, false) => {
                if r != q {
                    events.push(Event::Substitution { ref_pos, ref_base: r, read_base: q });
                }
                ref_pos += 1;
                i += 1;
            }
        }
    }
    events
}

/// Whether a substitution matches a base editor's chemistry.
fn accept_substitution(class: EventClass, ref_base: u8, read_base: u8) -> bool {
    let rb = ref_base.to_ascii_uppercase();
    let qb = read_base.to_ascii_uppercase();
    let is_abe = (rb == b'A' && qb == b'G') || (rb == b'T' && qb == b'C');
    let is_cbe = (rb == b'C' && qb == b'T') || (rb == b'G' && qb == b'A');
    match class {
        EventClass::AdenineBaseEdit => is_abe,
        EventClass::CytosineBaseEdit => is_cbe,
        EventClass::DualBaseEdit => is_abe || is_cbe,
        EventClass::DoubleStrandBreak | EventClass::PrimeEdit | EventClass::None => false,
    }
}

/// Call the events attributable to a single target from the already-extracted
/// alignment events. Returns the encoded per-target string, or `NONE` if
/// nothing overlaps its editing window.
pub fn call_target(all_events: &[Event], target_start: usize, target_type: &TargetType) -> String {
    let class = target_type.event_class();
    if class == EventClass::None || class == EventClass::PrimeEdit {
        return NONE_EVENT.to_string();
    }

    let (w0, w1) = target_type.editing_window();
    let win_start = target_start + w0;
    let win_end = target_start + w1;

    let kept: Vec<String> = all_events
        .iter()
        .filter(|ev| {
            let (es, ee) = ev.ref_span();
            let overlaps = es <= win_end && ee >= win_start;
            if !overlaps {
                return false;
            }
            match ev {
                Event::Deletion { .. } | Event::Insertion { .. } => {
                    class == EventClass::DoubleStrandBreak
                }
                Event::Substitution { ref_base, read_base, .. } => {
                    accept_substitution(class, *ref_base, *read_base)
                }
            }
        })
        .map(|ev| ev.encode())
        .collect();

    if kept.is_empty() {
        NONE_EVENT.to_string()
    } else {
        kept.join("&")
    }
}

/// Prime-edit haplotype classification written to the `pe` BAM tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimeEditCall {
    WildType,
    Precise,
    Partial,
    PrecisePlusByproduct,
    ScaffoldIncorporation,
    Indel,
    Other,
    NoCall,
}

impl PrimeEditCall {
    pub fn encode(&self) -> &'static str {
        match self {
            PrimeEditCall::WildType => "WT",
            PrimeEditCall::Precise => "PRECISE",
            PrimeEditCall::Partial => "PARTIAL",
            PrimeEditCall::PrecisePlusByproduct => "PRECISE_PLUS_BYPRODUCT",
            PrimeEditCall::ScaffoldIncorporation => "SCAFFOLD_INCORPORATION",
            PrimeEditCall::Indel => "INDEL",
            PrimeEditCall::Other => "OTHER",
            PrimeEditCall::NoCall => "NO_CALL",
        }
    }
}

/// Event strings plus optional target-aligned prime-edit classifications.
pub struct ReadEventCalls {
    pub events: String,
    pub prime_edits: Option<String>,
}

fn uppercase(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().map(|base| base.to_ascii_uppercase()).collect()
}

fn read_base_at_reference_position(
    reference_aligned: &[u8],
    read_aligned: &[u8],
    wanted_position: usize,
) -> Option<u8> {
    let mut ref_pos = 0usize;
    for (reference_base, read_base) in reference_aligned.iter().zip(read_aligned) {
        if *reference_base == GAP {
            continue;
        }
        if ref_pos == wanted_position {
            return if *read_base == GAP { None } else { Some(*read_base) };
        }
        ref_pos += 1;
    }
    None
}

/// Reconstruct the read haplotype spanning `[start, end)` in reference
/// coordinates. For a zero-width interval, capture insertions at `start`.
fn observed_haplotype(
    reference_aligned: &[u8],
    read_aligned: &[u8],
    start: usize,
    end: usize,
) -> Vec<u8> {
    let mut observed = Vec::new();
    let mut ref_pos = 0usize;

    for (reference_base, read_base) in reference_aligned.iter().zip(read_aligned) {
        if *reference_base == GAP {
            let insertion_is_in_range = if start == end {
                ref_pos == start
            } else {
                ref_pos >= start && ref_pos < end
            };
            if insertion_is_in_range && *read_base != GAP {
                observed.push(read_base.to_ascii_uppercase());
            }
            continue;
        }

        if ref_pos >= start && ref_pos < end && *read_base != GAP {
            observed.push(read_base.to_ascii_uppercase());
        }
        ref_pos += 1;
    }
    observed
}

fn prime_coordinates(
    reference_length: usize,
    target_start: usize,
    spec: &PrimeEditSpec,
) -> (usize, usize, usize, usize) {
    let edit_start = target_start + spec.edit_offset;
    let edit_end = edit_start + spec.reference.len();
    let window_start = edit_start.saturating_sub(spec.call_flank);
    let window_end = reference_length.min(edit_end.saturating_add(spec.call_flank));
    (edit_start, edit_end, window_start, window_end)
}

fn expected_haplotype(
    reference: &[u8],
    edit_start: usize,
    edit_end: usize,
    window_start: usize,
    window_end: usize,
    alternate: &[u8],
) -> Vec<u8> {
    let mut expected = uppercase(&reference[window_start..edit_start]);
    expected.extend(uppercase(alternate));
    expected.extend(uppercase(&reference[edit_end..window_end]));
    expected
}

fn partial_alleles(spec: &PrimeEditSpec) -> Vec<Vec<u8>> {
    let reference = uppercase(spec.reference.as_bytes());
    let alternate = uppercase(spec.alternate.as_bytes());
    let steps = reference.len().max(alternate.len());
    let mut partials = Vec::new();

    for incorporated in 1..steps {
        let candidate = match spec.strand {
            TargetStrand::Forward => {
                let mut candidate = alternate[..incorporated.min(alternate.len())].to_vec();
                candidate.extend_from_slice(&reference[incorporated.min(reference.len())..]);
                candidate
            }
            TargetStrand::Reverse => {
                let reference_cut = reference.len().saturating_sub(incorporated);
                let alternate_cut = alternate.len().saturating_sub(incorporated);
                let mut candidate = reference[..reference_cut].to_vec();
                candidate.extend_from_slice(&alternate[alternate_cut..]);
                candidate
            }
        };
        if candidate != reference && candidate != alternate && !partials.contains(&candidate) {
            partials.push(candidate);
        }
    }
    partials
}

fn contains_subsequence(haystack: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty() && haystack.windows(needle.len()).any(|window| window == needle)
}

/// Classify one prime-edit target against its programmed allele.
pub fn call_prime_edit(
    reference_aligned: &[u8],
    read_aligned: &[u8],
    reference_sequence: &[u8],
    target_start: usize,
    spec: &PrimeEditSpec,
) -> PrimeEditCall {
    let (edit_start, edit_end, window_start, window_end) =
        prime_coordinates(reference_sequence.len(), target_start, spec);

    let left_anchor = edit_start.checked_sub(1);
    let right_anchor = if edit_end < reference_sequence.len() {
        Some(edit_end)
    } else {
        None
    };
    let anchor_positions = [left_anchor, right_anchor];
    let anchors = anchor_positions.iter().filter_map(|anchor| *anchor);
    if anchors
        .map(|anchor| read_base_at_reference_position(reference_aligned, read_aligned, anchor))
        .any(|base| base.is_none())
    {
        return PrimeEditCall::NoCall;
    }

    let observed = observed_haplotype(
        reference_aligned,
        read_aligned,
        window_start,
        window_end,
    );
    let wild_type = uppercase(&reference_sequence[window_start..window_end]);
    let expected = expected_haplotype(
        reference_sequence,
        edit_start,
        edit_end,
        window_start,
        window_end,
        spec.alternate.as_bytes(),
    );

    if observed == wild_type {
        return PrimeEditCall::WildType;
    }
    if observed == expected {
        return PrimeEditCall::Precise;
    }
    if let Some(scaffold) = spec.scaffold_sequence.as_ref() {
        let scaffold = uppercase(scaffold.as_bytes());
        if contains_subsequence(&observed, &scaffold)
            && !contains_subsequence(&wild_type, &scaffold)
            && !contains_subsequence(&expected, &scaffold)
        {
            return PrimeEditCall::ScaffoldIncorporation;
        }
    }

    let observed_core = observed_haplotype(
        reference_aligned,
        read_aligned,
        edit_start,
        edit_end,
    );
    let alternate = uppercase(spec.alternate.as_bytes());
    if observed_core == alternate {
        return PrimeEditCall::PrecisePlusByproduct;
    }
    if partial_alleles(spec).contains(&observed_core) {
        return PrimeEditCall::Partial;
    }

    let has_indel = extract_all_events(reference_aligned, read_aligned)
        .iter()
        .any(|event| {
            let (start, end) = event.ref_span();
            start < window_end
                && end >= window_start
                && matches!(event, Event::Deletion { .. } | Event::Insertion { .. })
        });
    if has_indel {
        PrimeEditCall::Indel
    } else {
        PrimeEditCall::Other
    }
}

fn prime_event_string(
    all_events: &[Event],
    reference_length: usize,
    target_start: usize,
    spec: &PrimeEditSpec,
) -> String {
    let (_, _, window_start, window_end) =
        prime_coordinates(reference_length, target_start, spec);
    let events = all_events
        .iter()
        .filter(|event| {
            let (start, end) = event.ref_span();
            start < window_end && end >= window_start
        })
        .map(Event::encode)
        .collect::<Vec<_>>();
    if events.is_empty() {
        NONE_EVENT.to_string()
    } else {
        events.join("&")
    }
}

/// Call all target events for one aligned read against a reference record and
/// return the per-read event string: each target's call in target order,
/// joined by `_`. Returns an empty string when the reference declares no
/// targets.
pub fn call_read_event_details(
    reference_aligned: &[u8],
    read_aligned: &[u8],
    reference: &ReferenceRecord,
) -> ReadEventCalls {
    if reference.targets.is_empty() {
        return ReadEventCalls { events: String::new(), prime_edits: None };
    }

    let all = extract_all_events(reference_aligned, read_aligned);

    let locations = match reference.target_locations.as_ref() {
        Some(locs) => locs,
        // target_locations is filled during YAML load; if it somehow was not,
        // we cannot place the windows, so emit nothing rather than guess.
        None => return ReadEventCalls { events: String::new(), prime_edits: None },
    };

    let mut event_calls = Vec::with_capacity(reference.targets.len());
    let mut prime_calls = Vec::with_capacity(reference.targets.len());
    let mut has_prime_edit = false;

    for (target_index, (target_type, target_start)) in reference
        .target_types
        .iter()
        .zip(locations.iter())
        .enumerate()
    {
        if target_type == &TargetType::PrimeEdit {
            has_prime_edit = true;
            let spec = reference.prime_edits.get(&target_index).unwrap_or_else(|| {
                panic!("PrimeEdit target index {} has no specification", target_index)
            });
            event_calls.push(prime_event_string(
                &all,
                reference.sequence.len(),
                *target_start,
                spec,
            ));
            prime_calls.push(
                call_prime_edit(
                    reference_aligned,
                    read_aligned,
                    reference.sequence.as_bytes(),
                    *target_start,
                    spec,
                )
                .encode()
                .to_string(),
            );
        } else {
            event_calls.push(call_target(&all, *target_start, target_type));
            prime_calls.push("NA".to_string());
        }
    }

    ReadEventCalls {
        events: event_calls.join("_"),
        prime_edits: if has_prime_edit {
            Some(prime_calls.join("_"))
        } else {
            None
        },
    }
}

pub fn call_read_events(
    reference_aligned: &[u8],
    read_aligned: &[u8],
    reference: &ReferenceRecord,
) -> String {
    call_read_event_details(reference_aligned, read_aligned, reference).events
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_strategies::sequence_layout::ReferenceRecord;
    use std::collections::BTreeMap;

    fn prime_spec(reference: &str, alternate: &str, strand: TargetStrand) -> PrimeEditSpec {
        PrimeEditSpec {
            edit_offset: 6,
            reference: reference.to_string(),
            alternate: alternate.to_string(),
            strand,
            call_flank: 3,
            rtt_sequence: None,
            scaffold_sequence: None,
        }
    }

    fn alignment_with_replacement(
        reference: &[u8],
        edit_start: usize,
        reference_length: usize,
        observed_allele: &[u8],
    ) -> (Vec<u8>, Vec<u8>) {
        let mut reference_aligned = reference[..edit_start].to_vec();
        let mut read_aligned = reference[..edit_start].to_vec();
        let paired = reference_length.min(observed_allele.len());

        reference_aligned.extend_from_slice(&reference[edit_start..edit_start + paired]);
        read_aligned.extend_from_slice(&observed_allele[..paired]);
        if reference_length > paired {
            reference_aligned
                .extend_from_slice(&reference[edit_start + paired..edit_start + reference_length]);
            read_aligned.extend(std::iter::repeat(GAP).take(reference_length - paired));
        }
        if observed_allele.len() > paired {
            reference_aligned.extend(std::iter::repeat(GAP).take(observed_allele.len() - paired));
            read_aligned.extend_from_slice(&observed_allele[paired..]);
        }

        reference_aligned.extend_from_slice(&reference[edit_start + reference_length..]);
        read_aligned.extend_from_slice(&reference[edit_start + reference_length..]);
        (reference_aligned, read_aligned)
    }

    #[test]
    fn test_extract_perfect_match_is_empty() {
        let events = extract_all_events(b"ACGTACGT", b"ACGTACGT");
        assert!(events.is_empty());
    }

    #[test]
    fn test_extract_deletion() {
        //           A C G T A C G T
        // read:     A C G - - C G T   -> delete ref[3..=4], length 2 at pos 3
        let events = extract_all_events(b"ACGTACGT", b"ACG--CGT");
        assert_eq!(events, vec![Event::Deletion { length: 2, ref_pos: 3 }]);
        assert_eq!(events[0].encode(), "2D+3");
    }

    #[test]
    fn test_extract_insertion() {
        // ref:  A C - - G T   (ungapped ACGT)
        // read: A C T T G T   -> insert "TT" before ref pos 2
        let events = extract_all_events(b"AC--GT", b"ACTTGT");
        assert_eq!(events, vec![Event::Insertion { ref_pos: 2, bases: b"TT".to_vec() }]);
        assert_eq!(events[0].encode(), "2I+2+TT");
    }

    #[test]
    fn test_extract_substitution() {
        let events = extract_all_events(b"ACGT", b"AGGT");
        assert_eq!(
            events,
            vec![Event::Substitution { ref_pos: 1, ref_base: b'C', read_base: b'G' }]
        );
        assert_eq!(events[0].encode(), "1S+1+G");
    }

    #[test]
    fn test_extract_multiple_events() {
        // ref:  A C G T A C G T A C
        // read: A - G T A C G - - C  -> del at 1 (len1), del at 7 (len2)
        let events = extract_all_events(b"ACGTACGTAC", b"A-GTACG--C");
        assert_eq!(
            events,
            vec![
                Event::Deletion { length: 1, ref_pos: 1 },
                Event::Deletion { length: 2, ref_pos: 7 },
            ]
        );
    }

    #[test]
    fn test_dsb_keeps_indel_in_window_drops_substitution() {
        // Cas9WT window is (14,19). Build a deletion at pos 15 and a random
        // substitution at pos 16; only the deletion should survive.
        let events = vec![
            Event::Deletion { length: 3, ref_pos: 15 },
            Event::Substitution { ref_pos: 16, ref_base: b'A', read_base: b'C' },
        ];
        let call = call_target(&events, 0, &TargetType::Cas9WT);
        assert_eq!(call, "3D+15");
    }

    #[test]
    fn test_dsb_drops_indel_outside_window() {
        // Deletion at pos 2, well outside the Cas9WT (14,19) window.
        let events = vec![Event::Deletion { length: 3, ref_pos: 2 }];
        let call = call_target(&events, 0, &TargetType::Cas9WT);
        assert_eq!(call, NONE_EVENT);
    }

    #[test]
    fn test_abe_keeps_a_to_g_in_window_only() {
        // Cas9ABE window (2,19). A->G at pos 5 (kept), C->T at pos 6 (wrong
        // chemistry, dropped), A->G at pos 40 (out of window, dropped).
        let events = vec![
            Event::Substitution { ref_pos: 5, ref_base: b'A', read_base: b'G' },
            Event::Substitution { ref_pos: 6, ref_base: b'C', read_base: b'T' },
            Event::Substitution { ref_pos: 40, ref_base: b'A', read_base: b'G' },
        ];
        let call = call_target(&events, 0, &TargetType::Cas9ABE);
        assert_eq!(call, "1S+5+G");
    }

    #[test]
    fn test_abe_reverse_strand_t_to_c() {
        let events = vec![Event::Substitution { ref_pos: 10, ref_base: b'T', read_base: b'C' }];
        let call = call_target(&events, 0, &TargetType::Cas9ABE);
        assert_eq!(call, "1S+10+C");
    }

    #[test]
    fn test_cbe_keeps_c_to_t() {
        let events = vec![Event::Substitution { ref_pos: 8, ref_base: b'C', read_base: b'T' }];
        assert_eq!(call_target(&events, 0, &TargetType::Cas9CBE), "1S+8+T");
        // ABE should reject the same event
        assert_eq!(call_target(&events, 0, &TargetType::Cas9ABE), NONE_EVENT);
    }

    #[test]
    fn test_dual_base_edit_keeps_both() {
        let events = vec![
            Event::Substitution { ref_pos: 5, ref_base: b'A', read_base: b'G' },
            Event::Substitution { ref_pos: 9, ref_base: b'C', read_base: b'T' },
        ];
        let call = call_target(&events, 0, &TargetType::Cas9ABECBE);
        assert_eq!(call, "1S+5+G&1S+9+T");
    }

    #[test]
    fn test_base_editor_ignores_indels() {
        let events = vec![Event::Deletion { length: 2, ref_pos: 5 }];
        assert_eq!(call_target(&events, 0, &TargetType::Cas9ABE), NONE_EVENT);
    }

    #[test]
    fn test_static_target_is_none() {
        let events = vec![Event::Substitution { ref_pos: 1, ref_base: b'A', read_base: b'G' }];
        assert_eq!(call_target(&events, 0, &TargetType::Static), NONE_EVENT);
    }

    #[test]
    fn test_target_start_offsets_window() {
        // Same relative deletion, but target starts at reference position 100.
        // Cas9WT window (14,19) -> absolute [114,119]; deletion at 116 kept.
        let events = vec![Event::Deletion { length: 1, ref_pos: 116 }];
        assert_eq!(call_target(&events, 100, &TargetType::Cas9WT), "1D+116");
        // A deletion at 16 (would be in-window if target_start were 0) is now out.
        let events2 = vec![Event::Deletion { length: 1, ref_pos: 16 }];
        assert_eq!(call_target(&events2, 100, &TargetType::Cas9WT), NONE_EVENT);
    }

    fn make_reference(sequence: &str, targets: Vec<&str>, types: Vec<TargetType>) -> ReferenceRecord {
        let mut record = ReferenceRecord {
            sequence: sequence.to_string(),
            umi_configurations: BTreeMap::new(),
            targets: targets.iter().map(|t| t.to_string()).collect(),
            target_types: types,
            target_locations: None,
            prime_edits: BTreeMap::new(),
        };
        record.fill_and_validate_target_positions();
        record
    }

    #[test]
    fn test_call_read_events_end_to_end() {
        // Reference with one Cas9WT target. The target sequence sits at a known
        // offset, and we introduce a deletion inside its cut window.
        let refseq = "AAAAAAAAAAGGGGACGTACGTACGTACGTACGTAGGTTTTTTTTTT";
        // target is the 23bp protospacer starting at index 10:
        let target = &refseq[10..33];
        let reference = make_reference(refseq, vec![target], vec![TargetType::Cas9WT]);
        let target_start = reference.target_locations.as_ref().unwrap()[0];
        assert_eq!(target_start, 10);

        // Build a gapped alignment: perfect except a 2bp deletion at absolute
        // reference position 10+15 = 25 (inside window [24,29]).
        let ref_aligned: Vec<u8> = refseq.bytes().collect();
        let mut read_aligned = ref_aligned.clone();
        read_aligned[25] = GAP;
        read_aligned[26] = GAP;

        let call = call_read_events(&ref_aligned, &read_aligned, &reference);
        assert_eq!(call, "2D+25");
    }

    #[test]
    fn test_call_read_events_multiple_targets_joined() {
        // Two targets; only the second one is edited -> "NONE_<event>".
        let refseq = "AAAAAAAAAAGGGGACGTACGTACGTACGTACGTAGGTTTTTTTTTTCCCCGATCGATCGATCGATCGATCGGGGGGGGG";
        let t1 = &refseq[10..33];
        let t2 = &refseq[47..70];
        let reference = make_reference(
            refseq,
            vec![t1, t2],
            vec![TargetType::Cas9WT, TargetType::Cas9WT],
        );
        let t2_start = reference.target_locations.as_ref().unwrap()[1];

        let ref_aligned: Vec<u8> = refseq.bytes().collect();
        let mut read_aligned = ref_aligned.clone();
        // deletion inside t2's window [t2_start+14, t2_start+19]
        let del_pos = t2_start + 15;
        read_aligned[del_pos] = GAP;

        let call = call_read_events(&ref_aligned, &read_aligned, &reference);
        assert_eq!(call, format!("NONE_1D+{}", del_pos));
    }

    #[test]
    fn test_call_read_events_distinguishes_repeated_targets() {
        let target = "ACGTACGTACGTACGTACGTACG";
        let refseq = format!("{}TTTTTTTTTT{}", target, target);
        let reference = make_reference(
            &refseq,
            vec![target, target],
            vec![TargetType::Cas9WT, TargetType::Cas9WT],
        );
        assert_eq!(reference.target_locations, Some(vec![0, 33]));

        let ref_aligned = refseq.as_bytes().to_vec();
        let mut read_aligned = ref_aligned.clone();
        let deletion_position = 33 + 15;
        read_aligned[deletion_position] = GAP;

        assert_eq!(
            call_read_events(&ref_aligned, &read_aligned, &reference),
            format!("NONE_1D+{}", deletion_position)
        );
    }

    #[test]
    fn test_prime_edit_calls_wild_type_and_precise_substitution() {
        let reference = b"AAAACCCCGGGGTTTT";
        let spec = prime_spec("CC", "TT", TargetStrand::Forward);

        assert_eq!(
            call_prime_edit(reference, reference, reference, 0, &spec),
            PrimeEditCall::WildType
        );

        let (reference_aligned, read_aligned) =
            alignment_with_replacement(reference, 6, 2, b"TT");
        assert_eq!(
            call_prime_edit(&reference_aligned, &read_aligned, reference, 0, &spec),
            PrimeEditCall::Precise
        );
    }

    #[test]
    fn test_prime_edit_calls_precise_insertions_deletions_and_replacements() {
        let reference = b"AAAACCCCGGGGTTTT";
        for (reference_allele, alternate, observed) in [
            ("", "GG", b"GG".as_slice()),
            ("CC", "", b"".as_slice()),
            ("CC", "GGA", b"GGA".as_slice()),
        ] {
            let spec = prime_spec(reference_allele, alternate, TargetStrand::Forward);
            let (reference_aligned, read_aligned) = alignment_with_replacement(
                reference,
                6,
                reference_allele.len(),
                observed,
            );
            assert_eq!(
                call_prime_edit(&reference_aligned, &read_aligned, reference, 0, &spec),
                PrimeEditCall::Precise
            );
        }
    }

    #[test]
    fn test_prime_edit_partial_calls_follow_target_strand() {
        let reference = b"AAAACCCCGGGGTTTT";
        let forward = prime_spec("CC", "TT", TargetStrand::Forward);
        let reverse = prime_spec("CC", "TT", TargetStrand::Reverse);
        let (reference_aligned, forward_partial) =
            alignment_with_replacement(reference, 6, 2, b"TC");
        let (_, reverse_partial) = alignment_with_replacement(reference, 6, 2, b"CT");

        assert_eq!(
            call_prime_edit(&reference_aligned, &forward_partial, reference, 0, &forward),
            PrimeEditCall::Partial
        );
        assert_eq!(
            call_prime_edit(&reference_aligned, &reverse_partial, reference, 0, &reverse),
            PrimeEditCall::Partial
        );
    }

    #[test]
    fn test_prime_edit_calls_byproducts_scaffold_indels_and_no_call() {
        let reference = b"AAAACCCCGGGGTTTT";
        let mut spec = prime_spec("CC", "TT", TargetStrand::Forward);
        let (reference_aligned, mut precise_plus_byproduct) =
            alignment_with_replacement(reference, 6, 2, b"TT");
        precise_plus_byproduct[4] = b'T';
        assert_eq!(
            call_prime_edit(
                &reference_aligned,
                &precise_plus_byproduct,
                reference,
                0,
                &spec,
            ),
            PrimeEditCall::PrecisePlusByproduct
        );

        spec.scaffold_sequence = Some("AATT".to_string());
        let (scaffold_reference, scaffold_read) =
            alignment_with_replacement(reference, 6, 2, b"TTAATT");
        assert_eq!(
            call_prime_edit(&scaffold_reference, &scaffold_read, reference, 0, &spec),
            PrimeEditCall::ScaffoldIncorporation
        );

        let (indel_reference, indel_read) =
            alignment_with_replacement(reference, 6, 2, b"C");
        assert_eq!(
            call_prime_edit(&indel_reference, &indel_read, reference, 0, &spec),
            PrimeEditCall::Indel
        );

        let mut missing_anchor = reference.to_vec();
        missing_anchor[5] = GAP;
        assert_eq!(
            call_prime_edit(reference, &missing_anchor, reference, 0, &spec),
            PrimeEditCall::NoCall
        );
    }

    #[test]
    fn test_prime_edit_details_preserve_target_order_and_raw_events() {
        let reference_sequence = "AAAACCCCGGGGTTTT";
        let mut prime_edits = BTreeMap::new();
        prime_edits.insert(1, prime_spec("CC", "TT", TargetStrand::Forward));
        let reference = ReferenceRecord {
            sequence: reference_sequence.to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec!["AAAA".to_string(), reference_sequence.to_string()],
            target_types: vec![TargetType::Static, TargetType::PrimeEdit],
            target_locations: Some(vec![0, 0]),
            prime_edits,
        };
        let (reference_aligned, read_aligned) = alignment_with_replacement(
            reference_sequence.as_bytes(),
            6,
            2,
            b"TT",
        );

        let calls = call_read_event_details(&reference_aligned, &read_aligned, &reference);
        assert_eq!(calls.events, "NONE_1S+6+T&1S+7+T");
        assert_eq!(calls.prime_edits.as_deref(), Some("NA_PRECISE"));
    }

    #[test]
    fn test_call_read_events_no_targets_is_empty() {
        let reference = make_reference("ACGTACGTACGT", vec![], vec![]);
        assert_eq!(call_read_events(b"ACGTACGTACGT", b"ACGTACGTACGT", &reference), "");
    }
}
