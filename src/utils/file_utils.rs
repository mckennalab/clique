use flate2::bufread::GzDecoder;
use std::io::{BufRead, BufReader};
use std::fs::File;

pub fn get_reader(path: &str) -> Result<Box<dyn BufRead>, &'static str> {
    let file_type = path.split(".").collect::<Vec<&str>>().last().unwrap().clone();

    match file_type {
        "gz" => {
            let reader = Box::new(GzDecoder::new(BufReader::new(File::open(path).expect("Unable to open input known file"))));
            Ok(Box::new(BufReader::new(reader)))
        }
        _ => {
            Ok(Box::new(BufReader::new(File::open(path).expect("Unable to open known input file"))))
        }
    }
}
