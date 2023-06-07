use std::collections::{BTreeMap, HashMap, VecDeque};
use std::path::Path;
use crate::alignment_functions::align_reads;
use crate::InstanceLivedTempDir;
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration, UMISortType};
use crate::reference::fasta_reference::ReferenceManager;
use bio::io::fasta::*;
use itertools::Itertools;
use shardio::{Range, ShardReader, ShardWriter};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};
use crate::sequence_lookup::KnownLookup;
use rustc_hash::FxHashMap;
use crate::umis::sequence_clustering::{get_connected_components, input_list_to_graph, InputList, split_subgroup, string_distance};

pub fn collapse(reference: &String,
                final_output: &String,
                temp_directory: &InstanceLivedTempDir,
                read_structure: &SequenceLayoutDesign,
                max_reference_multiplier: &usize,
                min_read_length: &usize,
                read1: &String,
                read2: &String,
                index1: &String,
                index2: &String,
                threads: &usize,
                inversions: &bool) {


    // load up the reference files
    let rm = ReferenceManager::from(&reference, 8);

    // validate that each reference has the specified capture groups
    let validated_references = rm.references.iter().
        map(|rf| read_structure.validate_reference_sequence(&rf.sequence_u8)).all(|x| x == true);

    assert!(validated_references, "The reference sequences do not match the capture groups specified in the read structure file.");

    // align whatever combination of reads to the reference, collapsing down to a single sequence
    let aligned_temp = temp_directory.temp_file("aligned_reads.fasta");

    align_reads(&true,
                &false,
                &false,
                &rm,
                aligned_temp.as_path(),
                max_reference_multiplier,
                min_read_length,
                read1,
                read2,
                index1,
                index2,
                threads,
                inversions);

    // TODO: fix alignment to output as a sorted file and get rid of this rewrite step
    let mut sorted_input = output_fasta_to_sorted_shard_reader(&aligned_temp.as_path(), temp_directory, read_structure);

    // sort the reads by the tags
    read_structure.get_sorted_umi_configurations().iter().for_each(|tag| {
        match tag.sort_type {
            UMISortType::KnownTag => {
                sorted_input = sort_known_level(temp_directory, &sorted_input, &tag);
            },
            UMISortType::DegenerateTag => {
                sorted_input = sort_degenerate_level(temp_directory, &sorted_input, &tag);
            },
        }
    });

    // collapse the final reads down to a single sequence

    // write everything to the disk
}


fn consensus(input: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut consensus = Vec::new();

    for i in 0..input[0].len() {
        let mut counter = HashMap::new();

        input.iter().for_each(|vector| {
            *counter.entry(&vector[i]).or_insert(0) += 1;
        });

        let mut max = 0;
        let mut consensus_byte = b'N';

        for (byte, count) in counter {
            if count > max {
                max = count;
                consensus_byte = *byte;
            }
        }

        consensus.push(consensus_byte);
    }

    consensus
}

pub fn sort_degenerate_level(temp_directory: &InstanceLivedTempDir, reader: &ShardReader<SortingReadSetContainer>, tag: &UMIConfiguration) -> ShardReader<SortingReadSetContainer> {

    // first pass: read in the sequences we're clustering on, simply collecting the list. Then group this into clustered
    // tags that we 'correct' to the most parsimonious tag. Then, we can re-read the file and write out the corrected tags
    let reader2 = reader.clone();

    let mut input_list: FxHashMap<String, usize> = FxHashMap::default();

    reader2.iter_range(&Range::all()).unwrap().for_each(|x| {
        let mut sorting_read_set_container = x.unwrap();
        let next_key = sorting_read_set_container.ordered_unsorted_keys.pop_front().unwrap();
        assert_eq!(next_key.0, tag.symbol);
        *input_list.entry(FastaBase::to_string(&next_key.1)).or_insert(0) += 1;
    });

    let string_set = Vec::from_iter(input_list.keys()).iter().map(|s| s.as_bytes().to_vec()).collect::<Vec<Vec<u8>>>();
    let collection = InputList { strings: string_set, max_dist: u64::try_from(tag.max_distance).unwrap()};
    let graph = input_list_to_graph(&collection, string_distance, false);

    let cc = get_connected_components(&graph);
    println!("CC SIZE: {}", &cc.len());


    let mut final_correction: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

    for group in cc {
        let minilist = InputList { strings: group, max_dist: u64::try_from(tag.max_distance.clone()).unwrap()};
        let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

        let is_subgroups = split_subgroup(&mut minigraph);

        match is_subgroups {
            None => {
                let group = minilist.strings.clone();
                for s in minilist.strings {
                    final_correction.insert(s.clone(),consensus(&group));
                }
            }
            Some(x) => {
                for sgroup in &x {
                    for s in sgroup {
                        final_correction.insert(s.clone(),consensus(&sgroup));
                    }
                }
            }
        }
    }

    // create a new output
    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(aligned_temp, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();


    let mut sender = sharded_output.get_sender();

    reader.iter_range(&Range::all()).unwrap().for_each(|x| {
        let mut sorting_read_set_container = x.unwrap();
        let next_key = sorting_read_set_container.ordered_unsorted_keys.pop_front().unwrap();
        assert_eq!(next_key.0, tag.symbol);

        let corrected = final_correction.get(&FastaBase::to_vec_u8(&next_key.1)).unwrap();

        sorting_read_set_container.ordered_sorting_keys.push((next_key.0, FastaBase::from_vec_u8(corrected)));

        sender.send(sorting_read_set_container).unwrap();
    });


    sharded_output.finish().unwrap();

    let aligned_temp_ret = temp_directory.temp_file("first_pass_sorted.shareded");
    let reader = ShardReader::open(aligned_temp_ret).unwrap();

    reader
}

pub fn sort_known_level(temp_directory: &InstanceLivedTempDir, reader: &ShardReader<SortingReadSetContainer>, tag: &UMIConfiguration) -> ShardReader<SortingReadSetContainer> {

    // get the known lookup table
    let known_lookup = KnownLookup::from(tag);

    // create a new output
    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(aligned_temp, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();
    let mut sender = sharded_output.get_sender();

    reader.iter_range(&Range::all()).unwrap().for_each(|x| {
        let mut sorting_read_set_container = x.unwrap();
        let next_key = sorting_read_set_container.ordered_unsorted_keys.pop_front().unwrap();
        assert_eq!(next_key.0, tag.symbol);

        match known_lookup.correct(&FastaBase::to_string(&next_key.1), &tag.max_distance, false) {
            None => {}
            Some(x) => { sorting_read_set_container.ordered_sorting_keys.push((next_key.0, FastaBase::from_string(&x))) }
        }

        sender.send(sorting_read_set_container).unwrap();
    });


    sharded_output.finish().unwrap();

    let aligned_temp_ret = temp_directory.temp_file("first_pass_sorted.shareded");
    let reader = ShardReader::open(aligned_temp_ret).unwrap();

    reader
}

pub fn output_fasta_to_sorted_shard_reader(written_fasta_file: &Path,
                                           temp_directory: &InstanceLivedTempDir,
                                           read_structure: &SequenceLayoutDesign) -> ShardReader<SortingReadSetContainer> {
    let aligned_temp = temp_directory.temp_file("first_pass_sorted.sharded");

    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(aligned_temp, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();

    let mut sender = sharded_output.get_sender();

    let input_fasta = Reader::from_file(written_fasta_file).unwrap();

    let mut sorting_order = read_structure.umi_configurations.iter().map(|x| x.1.clone()).collect::<Vec<UMIConfiguration>>();
    sorting_order.sort_by(|a, b| a.order.cmp(&b.order));
    let sorted_tags = sorting_order.iter().map(|x| x.symbol).collect::<Vec<char>>();

    input_fasta.records().chunks(2).into_iter().for_each(|mut chunk| {
        let read_record = chunk.nth(1).unwrap().unwrap();
        let ref_record = chunk.nth(0).unwrap().unwrap();
        let read_tags = extract_output_tags(&read_record.id().to_string());
        let read_tags_ordered = VecDeque::from(sorted_tags.iter().
            map(|x| (x.clone(), read_tags.get(x).unwrap().as_bytes().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>())).collect::<Vec<(char, Vec<FastaBase>)>>());

        let new_sorted_read_container = SortingReadSetContainer {
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: read_tags_ordered,

            aligned_read: SortedAlignment {
                aligned_read: read_record.seq().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>(),
                aligned_ref: ref_record.seq().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>(),
                ref_name: "".to_string(),
            },
        };
        sender.send(new_sorted_read_container).unwrap();
    });

    sharded_output.finish().unwrap();

    let aligned_temp_ret = temp_directory.temp_file("first_pass_sorted.shareded");
    let reader = ShardReader::open(aligned_temp_ret).unwrap();

    reader
}

pub fn extract_output_tags(header_string: &String) -> BTreeMap<char, String> {
    assert!(header_string.starts_with(">"), "The header string does not start with a > character");
    let tokens = header_string.split("_").collect::<Vec<&str>>();
    assert!(tokens.len() >= 2, "The header string does not contain any tags");

    let mut tags: BTreeMap<char, String> = BTreeMap::new();
    tokens.get(1).unwrap().split(";").for_each(|tag| {
        let tag_tokens = tag.split("=").last().unwrap().split(":").collect::<Vec<&str>>();
        assert!(tag_tokens.len() == 2, "The tag {} does not contain a key and value", tag);
        let key = char::from(tag_tokens.get(0).unwrap().as_bytes()[0].clone());
        let value = tag_tokens.get(1).unwrap().to_string();
        println!("{} -> {}", key, value);
        tags.insert(key, value);
    });
    tags
}

//pub fn order_output_keys(header_string: &String, )

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use super::*;
    use crate::read_strategies::read_set::ReadIterator;

    #[test]
    fn test_extract_output_tags() {
        let test_string = String::from(">VH00708:168:AACK37GM5:1:1101:40992:1000_key=$:TCTCACGAGGTGGCTG;key=%:CGGTTCCGAAGT;key=^:AGGGTCTCGGCC");
        let tags = extract_output_tags(&test_string);
        assert_eq!(tags.len(), 3);
    }

    #[test]
    fn test_consensus() {
        let basic_seqs : Vec<Vec<u8>> = vec![String::from("ATCG").as_bytes().to_vec(), String::from("GCTA").as_bytes().to_vec(), String::from("ATCG").as_bytes().to_vec()];
        let consensus = consensus(&basic_seqs);
        assert_eq!(consensus, String::from("ATCG").as_bytes().to_vec());
    }


}