use std::fs::{remove_file, File};
use std::io::{prelude::*, BufReader, BufWriter, Lines};
use std::mem;

const BUFFER_CAPACITY: usize = 4_000_000_000;
const MAX_MEM_USE: usize = 4_000_000_000;

/*
fn main() {
    let filename = "logistic-chaos.txt";
    create_large_file(8_000_000_000, filename);
    extern_sort(filename, MAX_MEM_USE);
}
*/

pub(crate) fn create_tmp_list(file_count: usize) -> Vec<File> {
    (0..file_count).map(|c| tempfile::tempfile().unwrap()).into_iter().collect()
}


fn extern_sort(filename: &str, max_mem_use: usize) {
    let file = BufReader::with_capacity(BUFFER_CAPACITY, File::open(filename).unwrap());
    let mut v = vec![];
    let mut tmp_file_names = vec![];
    for x in file.lines() {
        v.push(x.unwrap().parse::<f64>().unwrap());
        if mem::size_of::<f64>() * v.len() > max_mem_use {
            sort_and_write_to_file(&mut v, &mut tmp_file_names);
        }
    }
    if v.len() > 0 {
        sort_and_write_to_file(&mut v, &mut tmp_file_names);
    }
    merge(&tmp_file_names, filename);
    clean_up(&tmp_file_names);
}

fn sort_and_write_to_file(v: &mut Vec<f64>, tmp_file_names: &mut Vec<String>) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let tmp_file_name = format!("tmp_sort_{}.txt", tmp_file_names.len());
    tmp_file_names.push(tmp_file_name.clone());
    println!("creating tmp file: {}", tmp_file_name);
    write_to_file(tmp_file_name, &v);
    v.clear();
}

fn clean_up(file_names: &Vec<String>) {
    file_names.iter().for_each(|p| remove_file(p).unwrap());
}

fn merge(tmp_file_names: &Vec<String>, file_name: &str) {
    println!("merging result ...");
    let result_file_name = format!("{}-sorted.txt", file_name.strip_suffix(".txt").unwrap());
    let mut file = BufWriter::with_capacity(BUFFER_CAPACITY, File::create(result_file_name).unwrap());
    let mut active_readers = tmp_file_names
        .iter()
        .map(|name| BufReader::with_capacity(BUFFER_CAPACITY, File::open(name).unwrap()).lines())
        .collect::<Vec<Lines<BufReader<File>>>>();
    let mut values = active_readers
        .iter_mut()
        .map(|r| r.next().unwrap().unwrap().parse::<f64>().unwrap())
        .collect::<Vec<f64>>();
    while active_readers.len() > 0 {
        let (i, max_val) = values
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        writeln!(file, "{}", max_val);
        if let Some(x) = active_readers[i].next() {
            values[i] = x.unwrap().parse::<f64>().unwrap();
        } else {
            values.remove(i);
            active_readers.remove(i);
        }
    }
    file.flush().unwrap();
}

fn write_to_file(filename: String, v: &Vec<f64>) {
    let mut buffer = BufWriter::with_capacity(BUFFER_CAPACITY, File::create(&filename).unwrap());
    for x in v.iter() {
        writeln!(buffer, "{}", x);
    }
    buffer.flush().unwrap();
}

fn create_large_file(size: usize, filename: &str) {
    println!("creating large file ...");
    let mut file = BufWriter::with_capacity(BUFFER_CAPACITY, File::create(filename).unwrap());
    let mut logistic = Logistic(0.35);
    let mut last_log = 0;
    let mut file_size = 0;
    while file_size < size {
        writeln!(file, "{}", logistic.next().unwrap());
        file_size += mem::size_of::<f64>();
        if file_size - last_log > 1_000_000 {
            println!("{}mb", file_size as f64 / 1_000_000.0);
            last_log = file_size;
        }
    }
    file.flush().unwrap();
}

struct Logistic(f64);
impl Iterator for Logistic {
    type Item = f64;
    fn next(&mut self) -> Option<Self::Item> {
        self.0 = 3.7 * self.0 * (1.0 - self.0);
        Some(self.0)
    }
}