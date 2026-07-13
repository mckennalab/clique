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

use crate::read_strategies::sequence_layout::{ReferenceRecord, TargetType};

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
        EventClass::DoubleStrandBreak | EventClass::None => false,
    }
}

/// Call the events attributable to a single target from the already-extracted
/// alignment events. Returns the encoded per-target string, or `NONE` if
/// nothing overlaps its editing window.
pub fn call_target(all_events: &[Event], target_start: usize, target_type: &TargetType) -> String {
    let class = target_type.event_class();
    if class == EventClass::None {
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

/// Call all target events for one aligned read against a reference record and
/// return the per-read event string: each target's call in target order,
/// joined by `_`. Returns an empty string when the reference declares no
/// targets.
pub fn call_read_events(
    reference_aligned: &[u8],
    read_aligned: &[u8],
    reference: &ReferenceRecord,
) -> String {
    if reference.targets.is_empty() {
        return String::new();
    }

    let all = extract_all_events(reference_aligned, read_aligned);

    let locations = match reference.target_locations.as_ref() {
        Some(locs) => locs,
        // target_locations is filled during YAML load; if it somehow was not,
        // we cannot place the windows, so emit nothing rather than guess.
        None => return String::new(),
    };

    reference
        .target_types
        .iter()
        .zip(locations.iter())
        .map(|(target_type, loc)| call_target(&all, *loc, target_type))
        .collect::<Vec<_>>()
        .join("_")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_strategies::sequence_layout::ReferenceRecord;
    use std::collections::BTreeMap;

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
    fn test_call_read_events_no_targets_is_empty() {
        let reference = make_reference("ACGTACGTACGT", vec![], vec![]);
        assert_eq!(call_read_events(b"ACGTACGTACGT", b"ACGTACGTACGT", &reference), "");
    }
}
