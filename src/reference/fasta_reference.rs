use fastq::OwnedRecord;
use bio::io::fasta::*;

pub struct Reference {
    pub sequence: Vec<u8>,
    pub name: Vec<u8>,
}

pub fn reference_file_to_struct(reference_file: &String) -> Reference {

    // reference loading / checking
    let mut reader = Reader::from_file(reference_file).unwrap();
    let fasta_entries: Vec<Record> = reader.records().map(|f| f.unwrap()).collect();
    assert_eq!(fasta_entries.len(), 1, "We can only run with single entry FASTA files");

    Reference {
        sequence: fasta_entries[0].seq().to_vec(),
        name: str::as_bytes(fasta_entries[0].id()).to_vec(),
    }
}