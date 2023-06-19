use mimalloc::MiMalloc;
use rust_htslib::{bam, bam::Read, bam::record::Aux};
use polars::{df, prelude::*};
use pyo3::prelude::*;
use pyo3_polars::{error::PyPolarsErr, PyDataFrame};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

fn main() {
    let bam_file = "/Volumes/archive/scratch/deardenlab/Jgilligan/Polistes_chi_meth/pass/bam_runid_b50030fdcc8d8f5caae75d8e3f6fd5fd011d9223_9_0.bam";
    let mut bam = bam::Reader::from_path(bam_file).unwrap();
    let header = bam::Header::from_template(bam.header());

    // Print out mod bases and probabilities using ML and MM tags
    for r in bam.records() {
        let record = r.unwrap();
        let seq = record.seq();
        let qual = record.qual();

        // let record = record.aux(b"ML");
        // println!("{:#?}", record); // ArrayU8

        let mm_record = record.aux(b"MM");
        if mm_record.is_err() {
            continue;
        }

        let mm_record = mm_record.unwrap();
        let mm_record = if let Aux::String(mm_record) = mm_record {
            mm_record
        } else {
            panic!("MM record is not a string");
        };

        // Drop semi-colon at the end
        let mm_record = &mm_record[..mm_record.len() - 1];

        // Collect into vec
        let mm_record: Vec<&str> = mm_record.split(",").collect();

        let modified_base_info = mm_record[0];
        let modified_base = modified_base_info.chars().nth(0).unwrap();

        // Drop first part
        let mm_record = &mm_record[1..];

        // Let's grab the ML Tag
        let ml_record = record.aux(b"ML");
        if ml_record.is_err() {
            panic!("MM Tag present but ML Tag is missing");
        }

        let ml_record = ml_record.unwrap();
        let ml_record = if let Aux::ArrayU8(ml_record) = ml_record {
            ml_record
        } else {
            panic!("ML record is not an ArrayU8");
        };

        // Convert to vec of u8
        let ml_record: Vec<u8> = ml_record.iter().map(|x| x).collect();

        // MM Tag is offset of modified base
        let mm_record: Vec<u16> = mm_record.iter().map(|x| x.parse::<u16>().unwrap()).collect();

        // MM tag tells us how many to skip
        // This is then followed by a comma separated list of how many seq bases of the stated base type to
        // skip, stored as a delta to the last and starting with 0 as the first (or next) base, starting from the
        // uncomplemented 5’ end of the SEQ field.

        // Iterate through possible modified bases, and print out the location and probabality, and modification type
        let mut mm_intermediate_count = 0;
        let mut mm_pos = 0;
        let mut total_modified_bases = 0;
        let mut total_possible_modified_bases = 0;
        let mut triplets: Vec<(u8, u8, u8)> = Vec::new();
        let mut triplet_watch = None;
        let mut triplet: (u8, u8, u8) = (0, 0, 0);

        println!("Modified base we are searching for: {}", modified_base);
        println!("Length: {}", seq.len());
        
        for (pos,x) in seq.as_bytes().iter().enumerate() {

            if triplet_watch == Some(1) {
                triplet.1 = *x;
                triplet_watch = Some(triplet_watch.unwrap() + 1);
            }

            if triplet_watch == Some(2) {
                triplet.2 = *x;
                triplet_watch = None;
                triplets.push(triplet);
            }

            if *x as char == modified_base {
                total_possible_modified_bases += 1;

                if mm_pos == mm_record.len() {
                    break;
                }

                if mm_intermediate_count == mm_record[mm_pos] {
                    // println!("{} {} {}", *x as char, pos, ml_record[mm_pos]);
                    mm_pos += 1;
                    mm_intermediate_count = 0;
                    total_modified_bases += 1;
                    triplet_watch = Some(1);
                    triplet.0 = *x;
                } else {
                    mm_intermediate_count += 1;
                }
            }
        }

        // Print up Seq ID and total number of modified bases, and total number of possible modified bases
        println!("Seq ID: {}", std::str::from_utf8(record.qname()).unwrap());
        println!("Total number of modified bases: {}", total_modified_bases);
        println!("Total number of possible modified bases: {}", total_possible_modified_bases);

        // Count all triples and print up abundances
        let mut triple_counts: Vec<(u8, u8, u8, u32)> = Vec::new();
        for triplet in triplets {
            let mut found = false;
            for (i, triple_count) in triple_counts.iter_mut().enumerate() {
                if triple_count.0 == triplet.0 && triple_count.1 == triplet.1 && triple_count.2 == triplet.2 {
                    triple_count.3 += 1;
                    found = true;
                    break;
                }
            }

            if !found {
                triple_counts.push((triplet.0, triplet.1, triplet.2, 1));
            }
        }

        // Print up triple counts
        for triple_count in triple_counts {
            println!("{} {} {} {}", triple_count.0 as char, triple_count.1 as char, triple_count.2 as char, triple_count.3);
        }


        

        /*
        Ok(
          String(
            "C+m?,1,1,3,3,1,0,3,15,11,3,3,0,2,0,2,1,10,0,9,3,3,7,16,12,3,3,2,0,2,20,7,2,3,2,11,13,3,3,2,0,2,20,3,3,2,2,0,3,13,0,13,1,2,6,0,1,1,11,10,3,4,3,2,0,2,13,10,2,3,2,0,2;",
          ),
        )
        */

        // let record = record.aux(b"MM");
        // println!("{:#?}", record);
    }

}
