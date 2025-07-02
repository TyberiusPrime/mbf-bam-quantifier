use crate::io::open_file;
use anyhow::{bail, Context, Result};
use std::{
    collections::{HashMap, HashSet},
    str::FromStr,
};

use crate::categorical::Categorical;

#[derive(Debug)]
pub struct GTFEntrys {
    pub seqname: Categorical,
    pub start: Vec<u64>,
    pub end: Vec<u64>,
    pub strand: Vec<Strand>,
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

    pub fn filter(&mut self, keep: &[bool]) {
        if keep.len() != self.start.len() {
            panic!("keep vector must be the same length as the GTFEntrys");
        }
        if keep.len() == 0 {
            return;
        }
        let mut iter = keep.iter();
        self.seqname.values.retain(|_| *iter.next().unwrap());

        let mut iter = keep.iter();
        self.start.retain(|_| *iter.next().unwrap());

        let mut iter = keep.iter();
        self.end.retain(|_| *iter.next().unwrap());

        let mut iter = keep.iter();
        self.strand.retain(|_| *iter.next().unwrap());

        for values in self.cat_attributes.values_mut() {
            let mut iter = keep.iter();
            values.values.retain(|_| *iter.next().unwrap());
        }
        for values in self.vec_attributes.values_mut() {
            let mut iter = keep.iter();
            assert_eq!(
                values.len(),
                keep.len(),
                "Vec attributes had different length than the rest of the GTFEntries:"
            );
            values.retain(|_| *iter.next().unwrap());
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

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Unstranded,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Unstranded => write!(f, "."),
        }
    }
}

impl FromStr for Strand {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            "." => Strand::Unstranded,
            "_" => Strand::Unstranded,
            _ => bail!("Invalid strand value: {}", s),
        })
    }
}

impl Into<i8> for &Strand {
    fn into(self) -> i8 {
        match self {
            Strand::Plus => 1,
            Strand::Minus => -1,
            Strand::Unstranded => 0,
        }
    }
}

//A gtf 'parser' (more like 'extract what we need') that doesn't
// allocate every line at least once
pub fn parse_minimal(
    filename: &str,
    accepted_features: HashSet<String>,
    accepted_tags: HashSet<String>,
) -> Result<HashMap<String, GTFEntrys>> {
    use linereader::LineReader; // non allocateding.
    let file = open_file(filename)?;
    let mut reader = LineReader::with_capacity(128 * 1024, file);
    let mut result = HashMap::new();

    while let Some(line) = reader.next_line() {
        let line = line.context("line reading error")?;
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b"#") {
            continue; // skip comments
        }
        let line = std::str::from_utf8(line).context("line was not utf8")?;
        if !line.ends_with("\n") {
            bail!("Line length exceed our buffer size, please increase the buffer size in the LineReader::with_capacity() call.");
        }
        let mut fields = line.split('\t');
        let seqname = fields.next().context("No seqname")?;
        let _source = fields.next().context("No source")?;
        let feature_type = fields.next().context("No feature type")?;
        if !accepted_features.contains(feature_type) {
            continue;
        }
        let start: u64 = fields
            .next()
            .context("No start")?
            .parse()
            .context("start not int")?;
        let start = start.checked_sub(1).context("start must be >= 1")?; // GTF is 1-based, convert to 0-based
        let end: u64 = fields
            .next()
            .context("No end")?
            .parse()
            .context("end not int")?;
        let _score = fields.next().context("no score")?;
        let strand: Strand = fields
            .next()
            .context("no strand")?
            .parse()
            .context("failed to parse strand. Allowed +-._")?;
        let _frame = fields.next().context("no frame")?;
        let attributes_str = fields.next().context("No attributes")?;
        let it = attributes_str
            .split_terminator(';')
            .map(str::trim_start)
            .filter(|x| !x.is_empty());
        let mut tags = Vec::new();
        for attr_value in it {
            let mut kv = attr_value.splitn(2, ' ');
            let key: &str = kv.next().unwrap();
            if !accepted_tags.contains(key) {
                continue;
            }
            let value: &str = kv.next().unwrap().trim_end().trim_matches('"');
            tags.push((key, value));
        }
        if !tags.is_empty() {
            let entry = result
                .entry(feature_type.to_string())
                .or_insert_with(GTFEntrys::new);
            entry.seqname.push(seqname);
            entry.start.push(start);
            entry.end.push(end);
            entry.strand.push(strand);
            let mut seen = HashSet::new();
            for (key, value) in tags {
                if seen.contains(key) {
                    bail!("doublicate attribute in GTF: {} in line: {}", key, line);
                }
                seen.insert(key);
                match entry.vec_attributes.entry(key.to_string()) {
                    std::collections::hash_map::Entry::Occupied(mut e) => {
                        e.get_mut().push(value.to_string());
                    }
                    std::collections::hash_map::Entry::Vacant(e) => {
                        e.insert(vector_new_empty_push(entry.count, value.to_string()));
                    }
                }
            }
            entry.count += 1;
        }
    }
    Ok(result)
}

/* pub fn parse_noodles_gtf(
    filename: &str,
    accepted_features: HashSet<String>,
    accepted_tags: HashSet<String>,
) -> Result<HashMap<String, GTFEntrys>> {
    use noodles_gff::feature::record::Attributes;
    let f = BufReader::new(open_file(filename)?);
    let mut reader = noodles_gtf::io::Reader::new(f);
    let mut out: HashMap<String, GTFEntrys> = HashMap::new();

    for result in reader.record_bufs() {
        let record = result?;

        let feature = std::str::from_utf8(record.ty()).unwrap();
        if !out.contains_key(feature) {
            if (!accepted_features.is_empty()) && (!accepted_features.contains(feature)) {
                continue;
            }
            let hm: GTFEntrys = GTFEntrys::new();
            out.insert(feature.to_string(), hm);
        }
        let target = out.get_mut(feature).unwrap();
        target
            .seqname
            .push(std::str::from_utf8(record.reference_sequence_name()).unwrap());
        target.start.push(record.start().get() as u64 - 1);
        target.end.push(record.end().get() as u64);
        target.strand.push(match record.strand() {
            noodles_gff::feature::record::Strand::Forward => Strand::Plus,
            noodles_gff::feature::record::Strand::Reverse => Strand::Minus,
            _ => Strand::Unstranded,
        });

        for res in record.attributes().iter() {
            let (key, value) = res?;
            let key = std::str::from_utf8(&key).context("Failed to parse key")?;
            if !accepted_tags.contains(key) {
                continue;
            }
            /* if (key.starts_with("gene") & (key != "gene_id") & (feature != "gene"))
                | (key.starts_with("transcript")
                    & (key != "transcript_id")
                    & (feature != "transcript"))
            {
                continue;
            } */
            let str_value = match &value {
                noodles_gff::feature::record::attributes::field::Value::String(value) => {
                    std::str::from_utf8(value).context("Failed to parse value")?
                }
                noodles_gff::feature::record::attributes::field::Value::Array(_) => continue,
            };

            if key.ends_with("_id") {
                // vec vs categorical seems to be almost performance neutral
                //just this push here (and I guess the fill-er-up below
                //takes about 3 seconds.
                target
                    .vec_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(str_value.to_string());
                    })
                    .unwrap_or_else(|| {
                        target.vec_attributes.insert(
                            key.to_string(),
                            vector_new_empty_push(target.count, str_value.to_string()),
                        );
                    });
            } else {
                // these tributes take about 1.5s to store (nd fill-er-up)
                target
                    .cat_attributes
                    .get_mut(key)
                    .map(|at| {
                        at.push(str_value);
                    })
                    .unwrap_or_else(|| {
                        target.cat_attributes.insert(
                            key.to_string(),
                            Categorical::new_empty_push(target.count, str_value),
                        );
                    });
            }
        }
    }

    Ok(out)
} */

/* pub fn parse_ensembl_gtf(
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
        let strand: Strand = strand.parse().context("Failet do parse strand")?;
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
} */
