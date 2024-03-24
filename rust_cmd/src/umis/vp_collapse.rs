use std::cmp::{max, Ordering};
use vpsearch::MetricSpace;
use crate::umis::known_list::FastaString;

/*

pub fn collapse_barcodes(strings: &mut [FastaString], max_distance: u32) -> Vec<FastaString> {
    let mut ret = Vec::new();
    match strings.len() {
        x if x > 2 => {
            let vantage_point_index = strings.len() - 1;
            let rest = &mut strings[..vantage_point_index];

            let half_idx = rest.len() / 2;

            let (near, far) = rest.split_at_mut(half_idx);
            let near_collapsed = collapse_barcodes(near, max_distance);
            let far_collapsed = collapse_barcodes(far, max_distance);
            ret
        },
        x if x == 2 => {
            // are we close enough?
            if strings[0].distance(&strings[1], &()) <= max_distance {
                let new_count = strings[0].count + strings[1].count;
                if strings[0].count > strings[1].count {
                    strings[0].count = new_count;
                    vec![strings[0].clone()]
                } else {
                    strings[1].count = new_count;
                    vec![strings[1].clone()]
                }
            } else {
                vec![strings[0].clone(),strings[1].clone()]
            }

        }
        x if x == 1 => {
            vec![strings[0].clone()]
        }
        _ => {
            panic!("0->-inf number of nodes")
        }
    }

}

fn sort_fasta_stings(vantage_point: &FastaString, strings: &mut [FastaString]) {
    for i in 0..strings.len() {
        strings[i].distance = vantage_point.distance(&strings[i], &());
    }
    strings.sort_unstable_by(|a, b| if a.distance < b.distance {Ordering::Less} else {Ordering::Greater});
}*/


pub struct MergingNode {
    radius: u32,

}