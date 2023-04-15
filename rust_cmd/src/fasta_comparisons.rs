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
