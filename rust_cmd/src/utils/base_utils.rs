use crate::fasta_comparisons::DEGENERATEBASES;

#[allow(dead_code)]
pub fn edit_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> usize {
    assert_eq!(str1.len(), str2.len());

    let mut dist: usize = 0;
    for i in 0..str1.len() {
        if !((DEGENERATEBASES.get(&str1[i]).is_some() && DEGENERATEBASES.get(&str1[i]).unwrap().contains_key(&str2[i])) ||
            (DEGENERATEBASES.get(&str2[i]).is_some() && DEGENERATEBASES.get(&str2[i]).unwrap().contains_key(&str1[i]))) {
            dist += 1;
        }
    }
    dist
}

pub fn is_valid_fasta_base(b: &u8) -> bool {
    matches!(b.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' |
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' |
        b'B' | b'D' | b'H' | b'V' | b'N'
    )
}

#[allow(dead_code)]
pub fn simple_edit_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> usize {
    assert_eq!(str1.len(), str2.len());

    let mut dist: usize = 0;
    for i in 0..str1.len() {
        if str1[i] != str2[i] {dist += 1}
    }
    dist
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edit_distance_identical() {
        let s1 = vec![b'A', b'C', b'G', b'T'];
        let s2 = vec![b'A', b'C', b'G', b'T'];
        assert_eq!(edit_distance(&s1, &s2), 0);
    }

    #[test]
    fn test_edit_distance_all_different() {
        let s1 = vec![b'A', b'A', b'A', b'A'];
        let s2 = vec![b'T', b'T', b'T', b'T'];
        assert_eq!(edit_distance(&s1, &s2), 4);
    }

    #[test]
    fn test_edit_distance_single_mismatch() {
        let s1 = vec![b'A', b'C', b'G', b'T'];
        let s2 = vec![b'A', b'C', b'G', b'A'];
        assert_eq!(edit_distance(&s1, &s2), 1);
    }

    #[test]
    fn test_edit_distance_degenerate_bases() {
        // R = A or G, so A vs R should be distance 0
        let s1 = vec![b'A'];
        let s2 = vec![b'R'];
        assert_eq!(edit_distance(&s1, &s2), 0);

        // N matches everything
        let s1 = vec![b'N'];
        let s2 = vec![b'T'];
        assert_eq!(edit_distance(&s1, &s2), 0);

        // Y = C or T, so G vs Y should be distance 1
        let s1 = vec![b'G'];
        let s2 = vec![b'Y'];
        assert_eq!(edit_distance(&s1, &s2), 1);
    }

    #[test]
    fn test_edit_distance_case_insensitive() {
        let s1 = vec![b'a'];
        let s2 = vec![b'A'];
        assert_eq!(edit_distance(&s1, &s2), 0);

        let s1 = vec![b'a'];
        let s2 = vec![b'a'];
        assert_eq!(edit_distance(&s1, &s2), 0);
    }

    #[test]
    #[should_panic]
    fn test_edit_distance_different_lengths() {
        let s1 = vec![b'A', b'C'];
        let s2 = vec![b'A'];
        edit_distance(&s1, &s2);
    }

    #[test]
    fn test_edit_distance_empty() {
        let s1: Vec<u8> = vec![];
        let s2: Vec<u8> = vec![];
        assert_eq!(edit_distance(&s1, &s2), 0);
    }

    #[test]
    fn test_is_valid_fasta_base_standard() {
        assert!(is_valid_fasta_base(&b'A'));
        assert!(is_valid_fasta_base(&b'C'));
        assert!(is_valid_fasta_base(&b'G'));
        assert!(is_valid_fasta_base(&b'T'));
        assert!(is_valid_fasta_base(&b'U'));
        assert!(is_valid_fasta_base(&b'N'));
    }

    #[test]
    fn test_is_valid_fasta_base_lowercase() {
        assert!(is_valid_fasta_base(&b'a'));
        assert!(is_valid_fasta_base(&b'c'));
        assert!(is_valid_fasta_base(&b'g'));
        assert!(is_valid_fasta_base(&b't'));
        assert!(is_valid_fasta_base(&b'n'));
    }

    #[test]
    fn test_is_valid_fasta_base_degenerate() {
        assert!(is_valid_fasta_base(&b'R'));
        assert!(is_valid_fasta_base(&b'Y'));
        assert!(is_valid_fasta_base(&b'S'));
        assert!(is_valid_fasta_base(&b'W'));
        assert!(is_valid_fasta_base(&b'K'));
        assert!(is_valid_fasta_base(&b'M'));
        assert!(is_valid_fasta_base(&b'B'));
        assert!(is_valid_fasta_base(&b'D'));
        assert!(is_valid_fasta_base(&b'H'));
        assert!(is_valid_fasta_base(&b'V'));
    }

    #[test]
    fn test_is_valid_fasta_base_invalid() {
        assert!(!is_valid_fasta_base(&b'-'));
        assert!(!is_valid_fasta_base(&b'X'));
        assert!(!is_valid_fasta_base(&b'0'));
        assert!(!is_valid_fasta_base(&b' '));
        assert!(!is_valid_fasta_base(&b'*'));
    }

    #[test]
    fn test_simple_edit_distance_identical() {
        let s1 = vec![b'A', b'C', b'G', b'T'];
        let s2 = vec![b'A', b'C', b'G', b'T'];
        assert_eq!(simple_edit_distance(&s1, &s2), 0);
    }

    #[test]
    fn test_simple_edit_distance_all_different() {
        let s1 = vec![b'A', b'A', b'A', b'A'];
        let s2 = vec![b'T', b'T', b'T', b'T'];
        assert_eq!(simple_edit_distance(&s1, &s2), 4);
    }

    #[test]
    fn test_simple_edit_distance_case_sensitive() {
        // simple_edit_distance is case-sensitive unlike edit_distance
        let s1 = vec![b'a'];
        let s2 = vec![b'A'];
        assert_eq!(simple_edit_distance(&s1, &s2), 1);
    }

    #[test]
    #[should_panic]
    fn test_simple_edit_distance_different_lengths() {
        let s1 = vec![b'A', b'C'];
        let s2 = vec![b'A'];
        simple_edit_distance(&s1, &s2);
    }
}
