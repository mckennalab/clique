use crate::alignment::fasta_bit_encoding::FastaString;

pub fn merge_two_reads(read1: FastaString, read2: FastaString) -> FastaString {

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_comp() {
        let fasta_string = FastaString::from("ACGT");
        let rev_comp = fasta_string.reverse_complement();
        assert_eq!(fasta_string, rev_comp);
    }
}