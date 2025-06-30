use crate::io::open_file;
use anyhow::{Context, Result};
use std::{
    collections::{HashMap, HashSet},
    io::{BufRead, BufReader},
};

use crate::categorical::Categorical;

#[derive(Debug)]
pub struct GTFEntrys {
    pub seqname: Categorical,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub strand: Vec<i8>,
    pub cat_attributes: HashMap<String, Categorical>,
    pub vec_attributes: HashMap<String, Vec<String>>,
    count: u32,
}

impl GTFEntrys {
    pub fn new() -> GTFEntrys {
        GTFEntrys {
            seqname: Categorical::new(),
            start: Vec::new(),
            end: Vec::new(),
            strand: Vec::new(),
            cat_attributes: HashMap::new(),
            vec_attributes: HashMap::new(),
            count: 0,
        }
    }

    /* pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    } */
}

/// a helper that creates a vector, fills it with empty strings up to count
/// then adds value
/// similar to Categorical.new_empty_push
fn vector_new_empty_push(count: u32, value: String) -> Vec<String> {
    let mut res = Vec::new();
    res.resize(count as usize, "".to_string());
    res.push(value);
    res
}

pub fn parse_ensembl_gtf(
    filename: &str,
    accepted_features: HashSet<String>,
) -> Result<HashMap<String, GTFEntrys>> {
    // this is good but it still iterates through parts of the input
    // three times!
    let f = BufReader::new(open_file(filename)?);
    let mut out: HashMap<String, GTFEntrys> = HashMap::new();
    for line in f.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(9, '\t');
        let seqname = parts.next().context("Failed to find seqname")?;
        parts.next(); //consume source
        let feature = parts.next().context("Failed to find feature")?;
        if !out.contains_key(feature) {
            if (!accepted_features.is_empty()) && (!accepted_features.contains(feature)) {
                continue;
            }
            let hm: GTFEntrys = GTFEntrys::new();
            out.insert(feature.to_string(), hm);
        }
        let start: u64 = parts.next().context("Failed to find start")?.parse()?;
        let start = start - 1;
        let end: u64 = parts.next().context("Failed to find start")?.parse()?;
        parts.next(); //consume score
        let strand = parts.next().context("Failed to find start")?;
        let strand: i8 = if strand == "+" {
            1
        } else if strand == "-" {
            -1
        } else {
            0
        };
        let target = out.get_mut(feature).unwrap();
        target.seqname.push(seqname);
        target.start.push(start);
        target.end.push(end);
        target.strand.push(strand);
        let mut tag_count = 0;
        parts.next(); //consume frame
        let attributes = parts.next().context("Failed to find attributes")?;
        let it = attributes
            .split_terminator(';')
            .map(str::trim_start)
            .filter(|x| !x.is_empty());
        for attr_value in it {
            let mut kv = attr_value.splitn(2, ' ');
            let mut key: &str = kv.next().unwrap();
            if key == "tag" {
                if feature != "transcript" {
                    // only transcripts have tags!
                    continue;
                }
                if tag_count == 0 {
                    key = "tag0"
                } else if tag_count == 1 {
                    key = "tag1"
                } else if tag_count == 2 {
                    key = "tag2"
                } else if tag_count == 3 {
                    key = "tag3"
                } else if tag_count == 4 {
                    key = "tag4"
                } else if tag_count == 5 {
                    key = "tag5"
                } else {
                    continue; // silently swallow further tags
                }
                tag_count += 1;
            }
            if (key.starts_with("gene") & (key != "gene_id") & (feature != "gene"))
                | (key.starts_with("transcript")
                    & (key != "transcript_id")
                    & (feature != "transcript"))
            {
                continue;
            }
            let value: &str = kv.next().unwrap().trim_matches('"');
            if key.ends_with("_id") {
                // vec vs categorical seems to be almost performance neutral
                //just this push here (and I guess the fill-er-up below
                //takes about 3 seconds.
                target
                    .vec_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(value.to_string());
                    })
                    .unwrap_or_else(|| {
                        target.vec_attributes.insert(
                            key.to_string(),
                            vector_new_empty_push(target.count, value.to_string()),
                        );
                    });
            } else {
                // these tributes take about 1.5s to store (nd fill-er-up)
                target
                    .cat_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(value);
                    })
                    .unwrap_or_else(|| {
                        target.cat_attributes.insert(
                            key.to_string(),
                            Categorical::new_empty_push(target.count, value),
                        );
                    });
            }
        }
        target.count += 1;
        for (_key, value) in target.cat_attributes.iter_mut() {
            if (value.len() as u32) < target.count {
                value.push("");
            }
        }
        for (_key, value) in target.vec_attributes.iter_mut() {
            if (value.len() as u32) < target.count {
                value.push("".to_string());
            }
        }
    }

    Ok(out)
}
