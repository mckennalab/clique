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
