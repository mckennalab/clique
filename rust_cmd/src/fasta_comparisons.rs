use std::{collections::HashMap, hash::BuildHasherDefault};
use nohash_hasher::NoHashHasher;

// sets of known characters: our standard DNA alphabet and a second version with known gaps.
// These are used to mask known values when looking for extractable UMI/ID/barcode sequences.
// Also mappings from degenerate FASTA bases to their possible ACGT values.
lazy_static! {
    pub static ref KNOWNBASES: HashMap::<u8, u8, BuildHasherDefault<NoHashHasher<u8>>> = {
        let mut hashedvalues: HashMap::<u8, u8, BuildHasherDefault<NoHashHasher<u8>>> = HashMap::with_capacity_and_hasher(8, BuildHasherDefault::default());
        hashedvalues.insert(b'a', b'A');
        hashedvalues.insert(b'A', b'A');
        hashedvalues.insert(b'c', b'C');
        hashedvalues.insert(b'C', b'C');
        hashedvalues.insert(b'g', b'G');
        hashedvalues.insert(b'G', b'G');
        hashedvalues.insert(b't', b'T');
        hashedvalues.insert(b'T', b'T');
        hashedvalues
    };

    pub static ref DEGENERATEBASES: HashMap::<u8, HashMap<u8, bool>> = {
        let mut hashedvalues = HashMap::new(); // with_capacity_and_hasher(15, BuildHasherDefault::default());
        hashedvalues.insert(b'A', HashMap::from([('A' as u8, true), ('a' as u8, true)]));
        hashedvalues.insert(b'a', HashMap::from([('A' as u8, true), ('a' as u8, true)]));

        hashedvalues.insert(b'C', HashMap::from([('C' as u8, true), ('c' as u8, true)]));
        hashedvalues.insert(b'c', HashMap::from([('C' as u8, true), ('c' as u8, true)]));

        hashedvalues.insert(b'G', HashMap::from([('G' as u8, true), ('g' as u8, true)]));
        hashedvalues.insert(b'g', HashMap::from([('G' as u8, true), ('g' as u8, true)]));

        hashedvalues.insert(b'T', HashMap::from([('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b't', HashMap::from([('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'R', HashMap::from([('A' as u8, true), ('a' as u8, true), ('G' as u8, true), ('g' as u8, true)]));
        hashedvalues.insert(b'r', HashMap::from([('A' as u8, true), ('a' as u8, true), ('G' as u8, true), ('g' as u8, true)]));

        hashedvalues.insert(b'Y', HashMap::from([('C' as u8, true), ('c' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'y', HashMap::from([('C' as u8, true), ('c' as u8, true), ('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'K', HashMap::from([('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'k', HashMap::from([('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'M', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true)]));
        hashedvalues.insert(b'm', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true)]));

        hashedvalues.insert(b'S', HashMap::from([('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true)]));
        hashedvalues.insert(b's', HashMap::from([('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true)]));

        hashedvalues.insert(b'W', HashMap::from([('A' as u8, true), ('a' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'w', HashMap::from([('A' as u8, true), ('a' as u8, true), ('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'B', HashMap::from([('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'b', HashMap::from([('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'D', HashMap::from([('A' as u8, true), ('a' as u8, true), ('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'd', HashMap::from([('A' as u8, true), ('a' as u8, true), ('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'H', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'h', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true), ('T' as u8, true), ('t' as u8, true)]));

        hashedvalues.insert(b'V', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true)]));
        hashedvalues.insert(b'v', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true)]));

        hashedvalues.insert(b'N', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues.insert(b'n', HashMap::from([('A' as u8, true), ('a' as u8, true), ('C' as u8, true), ('c' as u8, true), ('G' as u8, true), ('g' as u8, true), ('T' as u8, true), ('t' as u8, true)]));
        hashedvalues
    };

    pub static ref KNOWNBASESPLUSINSERT: HashMap::<u8, u8, BuildHasherDefault<NoHashHasher<u8>>> = {
        let mut hashedvalues = HashMap::with_capacity_and_hasher(9, BuildHasherDefault::default());
        hashedvalues.insert(b'a', b'A');
        hashedvalues.insert(b'A', b'A');
        hashedvalues.insert(b'c', b'C');
        hashedvalues.insert(b'C', b'C');
        hashedvalues.insert(b'g', b'G');
        hashedvalues.insert(b'G', b'G');
        hashedvalues.insert(b't', b'T');
        hashedvalues.insert(b'T', b'T');
        hashedvalues.insert(b'-', b'-');
        hashedvalues
    };

    pub static ref REVERSECOMP: HashMap::<u8, u8, BuildHasherDefault<NoHashHasher<u8>>> = {
            let mut hashedvalues = HashMap::with_capacity_and_hasher(8, BuildHasherDefault::default());
            hashedvalues.insert(b'a', b'T');
            hashedvalues.insert(b'A', b'T');
            hashedvalues.insert(b'c', b'G');
            hashedvalues.insert(b'C', b'G');
            hashedvalues.insert(b'g', b'C');
            hashedvalues.insert(b'G', b'C');
            hashedvalues.insert(b't', b'A');
            hashedvalues.insert(b'T', b'A');
            hashedvalues
        };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_knownbases_standard() {
        assert_eq!(KNOWNBASES.get(&b'A'), Some(&b'A'));
        assert_eq!(KNOWNBASES.get(&b'a'), Some(&b'A'));
        assert_eq!(KNOWNBASES.get(&b'C'), Some(&b'C'));
        assert_eq!(KNOWNBASES.get(&b'c'), Some(&b'C'));
        assert_eq!(KNOWNBASES.get(&b'G'), Some(&b'G'));
        assert_eq!(KNOWNBASES.get(&b'g'), Some(&b'G'));
        assert_eq!(KNOWNBASES.get(&b'T'), Some(&b'T'));
        assert_eq!(KNOWNBASES.get(&b't'), Some(&b'T'));
    }

    #[test]
    fn test_knownbases_excludes_degenerate() {
        assert_eq!(KNOWNBASES.get(&b'N'), None);
        assert_eq!(KNOWNBASES.get(&b'R'), None);
        assert_eq!(KNOWNBASES.get(&b'-'), None);
    }

    #[test]
    fn test_knownbasesplusinsert_includes_gap() {
        assert_eq!(KNOWNBASESPLUSINSERT.get(&b'-'), Some(&b'-'));
        assert_eq!(KNOWNBASESPLUSINSERT.get(&b'A'), Some(&b'A'));
        assert_eq!(KNOWNBASESPLUSINSERT.get(&b'N'), None);
    }

    #[test]
    fn test_reversecomp_standard() {
        assert_eq!(REVERSECOMP.get(&b'A'), Some(&b'T'));
        assert_eq!(REVERSECOMP.get(&b'a'), Some(&b'T'));
        assert_eq!(REVERSECOMP.get(&b'T'), Some(&b'A'));
        assert_eq!(REVERSECOMP.get(&b't'), Some(&b'A'));
        assert_eq!(REVERSECOMP.get(&b'G'), Some(&b'C'));
        assert_eq!(REVERSECOMP.get(&b'g'), Some(&b'C'));
        assert_eq!(REVERSECOMP.get(&b'C'), Some(&b'G'));
        assert_eq!(REVERSECOMP.get(&b'c'), Some(&b'G'));
    }

    #[test]
    fn test_reversecomp_excludes_others() {
        assert_eq!(REVERSECOMP.get(&b'N'), None);
        assert_eq!(REVERSECOMP.get(&b'-'), None);
    }

    #[test]
    fn test_degeneratebases_standard_bases() {
        // A matches A/a
        assert!(DEGENERATEBASES.get(&b'A').unwrap().contains_key(&b'A'));
        assert!(DEGENERATEBASES.get(&b'A').unwrap().contains_key(&b'a'));
        assert!(!DEGENERATEBASES.get(&b'A').unwrap().contains_key(&b'C'));
    }

    #[test]
    fn test_degeneratebases_r_purine() {
        // R = A or G
        let r_map = DEGENERATEBASES.get(&b'R').unwrap();
        assert!(r_map.contains_key(&b'A'));
        assert!(r_map.contains_key(&b'a'));
        assert!(r_map.contains_key(&b'G'));
        assert!(r_map.contains_key(&b'g'));
        assert!(!r_map.contains_key(&b'C'));
        assert!(!r_map.contains_key(&b'T'));
    }

    #[test]
    fn test_degeneratebases_y_pyrimidine() {
        // Y = C or T
        let y_map = DEGENERATEBASES.get(&b'Y').unwrap();
        assert!(y_map.contains_key(&b'C'));
        assert!(y_map.contains_key(&b'T'));
        assert!(!y_map.contains_key(&b'A'));
        assert!(!y_map.contains_key(&b'G'));
    }

    #[test]
    fn test_degeneratebases_n_any() {
        // N = A, C, G, or T
        let n_map = DEGENERATEBASES.get(&b'N').unwrap();
        assert!(n_map.contains_key(&b'A'));
        assert!(n_map.contains_key(&b'C'));
        assert!(n_map.contains_key(&b'G'));
        assert!(n_map.contains_key(&b'T'));
        assert_eq!(n_map.len(), 8); // 4 bases * 2 cases each
    }

    #[test]
    fn test_degeneratebases_case_insensitive_keys() {
        // Lowercase keys should have same mappings as uppercase
        let upper = DEGENERATEBASES.get(&b'R').unwrap();
        let lower = DEGENERATEBASES.get(&b'r').unwrap();
        assert_eq!(upper.len(), lower.len());
        for key in upper.keys() {
            assert!(lower.contains_key(key));
        }
    }

    #[test]
    fn test_degeneratebases_all_iupac_codes_present() {
        // All IUPAC degenerate codes should be present (both cases)
        let codes = [b'A', b'C', b'G', b'T', b'R', b'Y', b'K', b'M', b'S', b'W', b'B', b'D', b'H', b'V', b'N'];
        for code in codes {
            assert!(DEGENERATEBASES.contains_key(&code), "Missing uppercase code {}", code as char);
            assert!(DEGENERATEBASES.contains_key(&code.to_ascii_lowercase()), "Missing lowercase code {}", code as char);
        }
    }

    #[test]
    fn test_degeneratebases_b_not_a() {
        // B = C, G, T (not A)
        let b_map = DEGENERATEBASES.get(&b'B').unwrap();
        assert!(b_map.contains_key(&b'C'));
        assert!(b_map.contains_key(&b'G'));
        assert!(b_map.contains_key(&b'T'));
        assert!(!b_map.contains_key(&b'A'));
    }

    #[test]
    fn test_degeneratebases_d_not_c() {
        // D = A, G, T (not C)
        let d_map = DEGENERATEBASES.get(&b'D').unwrap();
        assert!(d_map.contains_key(&b'A'));
        assert!(d_map.contains_key(&b'G'));
        assert!(d_map.contains_key(&b'T'));
        assert!(!d_map.contains_key(&b'C'));
    }

    #[test]
    fn test_degeneratebases_h_not_g() {
        // H = A, C, T (not G)
        let h_map = DEGENERATEBASES.get(&b'H').unwrap();
        assert!(h_map.contains_key(&b'A'));
        assert!(h_map.contains_key(&b'C'));
        assert!(h_map.contains_key(&b'T'));
        assert!(!h_map.contains_key(&b'G'));
    }

    #[test]
    fn test_degeneratebases_v_not_t() {
        // V = A, C, G (not T)
        let v_map = DEGENERATEBASES.get(&b'V').unwrap();
        assert!(v_map.contains_key(&b'A'));
        assert!(v_map.contains_key(&b'C'));
        assert!(v_map.contains_key(&b'G'));
        assert!(!v_map.contains_key(&b'T'));
    }
}
