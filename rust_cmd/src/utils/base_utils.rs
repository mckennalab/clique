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

#[allow(dead_code)]
pub fn simple_edit_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> usize {
    assert_eq!(str1.len(), str2.len());

    let mut dist: usize = 0;
    for i in 0..str1.len() {
        if str1[i] != str2[i] {dist += 1}
    }
    dist
}
