use anyhow::Result;
use super::chunked_genome::{Chunk, ChunkedGenome};
use super::{add_hashmaps, OurTree};
use crate::bam_ext::{open_bam, BamRecordExtensions};
use rust_htslib::bam::{record::Aux, Read};
use rayon::prelude::*;
use rust_htslib::bam;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, RwLock};


/// python wrapper for py_count_reads_unstranded
pub fn count_reads_unstranded(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
    each_read_counts_once: bool,
) -> Result<HashMap<String, u32>> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new(gene_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();
                let (tree, gene_ids) = trees.get(&chunk.chr).unwrap();

                let counts = count_reads_in_region_unstranded(
                    bam,
                    //&chunk.tree,
                    tree,
                    chunk.tid,
                    chunk.start,
                    chunk.stop,
                    gene_ids.len() as u32,
                    each_read_counts_once,
                );
                let mut total = 0;
                let mut outside = 0;
                let mut result: HashMap<String, u32> = match counts {
                    Ok(counts) => {
                        let mut res = HashMap::new();
                        for (gene_no, cnt) in counts.0.iter().enumerate() {
                            let gene_id = &gene_ids[gene_no];
                            res.insert(gene_id.to_string(), *cnt);
                            total += cnt;
                        }
                        outside += counts.1;
                        res
                    }
                    _ => HashMap::new(),
                };
                result.insert("_total".to_string(), total);
                result.insert("_outside".to_string(), outside);
                result.insert(format!("_{}", chunk.chr), total);
                result
            })
            .reduce(HashMap::<String, u32>::new, add_hashmaps);
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        Ok(result)
    })
}

fn to_hashmap(counts: Vec<u32>, gene_ids: &Vec<String>, chr: &str) -> HashMap<String, u32> {
    let mut total = 0;
    let mut result = HashMap::new();
    for (gene_no, cnt) in counts.iter().enumerate() {
        let gene_id = &gene_ids[gene_no];
        *result.entry(gene_id.to_string()).or_insert(0) += *cnt;
        total += cnt;
    }
    result.insert("_total".to_string(), total);
    result.insert(format!("_{}", chr), total);
    result
}

impl ReadCatcher for Arc<RwLock<HashSet<String>>> {
    fn catch(&mut self, read: &bam::Record) {
        self.write()
            .unwrap()
            .insert(std::str::from_utf8(read.qname()).unwrap().to_string());
    }
}

/// python wrapper for py_count_reads_stranded
pub fn count_reads_stranded(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
    gene_trees: HashMap<String, (OurTree, Vec<String>)>,
    each_read_counts_once: bool,
    matching_read_output_bam_filename: Option<&str>, //todo: remove this? or is this actually
    //useful for debugging?
) -> Result<(HashMap<String, u32>, HashMap<String, u32>)> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new(gene_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let matching_read_catcher = match &matching_read_output_bam_filename {
        None => None,
        Some(_) => Some(Arc::new(RwLock::new(HashSet::new()))),
    };

    let result = pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();
                let (tree, gene_ids) = trees.get(&chunk.chr).unwrap();

                let both_counts = match &matching_read_catcher {
                    None => count_reads_in_region_stranded(
                        bam,
                        tree,
                        chunk.tid,
                        chunk.start,
                        chunk.stop,
                        gene_ids.len() as u32,
                        each_read_counts_once,
                        (),
                    ),
                    Some(catcher) => count_reads_in_region_stranded(
                        bam,
                        tree,
                        chunk.tid,
                        chunk.start,
                        chunk.stop,
                        gene_ids.len() as u32,
                        each_read_counts_once,
                        catcher.clone(),
                    ),
                };

                let both_counts = both_counts.unwrap_or_else(|_| (Vec::new(), Vec::new(), 0));

                let mut result = (
                    to_hashmap(both_counts.0, gene_ids, &chunk.chr),
                    to_hashmap(both_counts.1, gene_ids, &chunk.chr),
                );
                result.0.insert("_outside".to_string(), both_counts.2);
                result
            })
            .reduce(
                || (HashMap::<String, u32>::new(), HashMap::<String, u32>::new()),
                add_dual_hashmaps,
            );
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        //TODO:
        //write actually to bam if catcher..
        Ok(result)
    });
    match matching_read_catcher {
        Some(catcher) => {
            let catcher = Arc::try_unwrap(catcher).unwrap().into_inner().unwrap();
            let output_filename = matching_read_output_bam_filename.unwrap();
            let mut bam_in = open_bam(filename, index_filename).unwrap();
            bam_in.fetch(bam::FetchDefinition::All)?;

            let header = bam::Header::from_template(bam_in.header());
            {
                let mut bam_out =
                    bam::Writer::from_path(output_filename, &header, bam::Format::Bam)?;
                for read in bam_in.rc_records() {
                    let read = read?; //if it fails her something is *very* wrong
                    if catcher.contains(std::str::from_utf8(read.qname()).unwrap()) {
                        bam_out.write(&read)?;
                    }
                }
            }
            bam::index::build(output_filename, None, bam::index::Type::Bai, 4)
                .expect("Failed to build bam index");
        }
        None => {}
    };
    result
}

fn add_dual_hashmaps(
    a: (HashMap<String, u32>, HashMap<String, u32>),
    b: (HashMap<String, u32>, HashMap<String, u32>),
) -> (HashMap<String, u32>, HashMap<String, u32>) {
    (add_hashmaps(a.0, b.0), add_hashmaps(a.1, b.1))
}

/// count_reads_in_region_unstranded
///
/// counts the unstranded reads in a region,
/// matching the to the tree entries.
///
/// if each_read_counts_once is set, each read can count only for one gene - the first
/// one where a block hits
/// otherwise reads having blocks that can not be assigned to just one gene count for
/// all they're hitting
fn count_reads_in_region_unstranded(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    gene_count: u32,
    each_read_counts_once: bool,
) -> Result<(Vec<u32>, u32)> {
    let mut result = vec![0; gene_count as usize];
    let mut multimapper_dedup: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut gene_nos_seen = HashSet::<u32>::new();
    let mut outside_count = 0;
    let mut read: bam::Record = bam::Record::new();
    bam.fetch((tid, start as u64, stop as u64))?;
    while let Some(bam_result) = bam.read(&mut read) {
        bam_result?;
        // do not count multiple blocks matching in one gene multiple times
        gene_nos_seen.clear();
        let mut hit = false;
        let mut skipped = false;
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            skipped = true;
        }
        if !skipped {
            let blocks = read.blocks();
            for iv in blocks.iter() {
                if (iv.1 < start) || iv.0 >= stop || ((iv.0 < start) && (iv.1 >= start)) {
                    // if this block is outside of the region
                    // don't count it at all.
                    // if it is on a block boundary
                    // only count it for the left side.
                    // which is ok, since we place the blocks to the right
                    // of our intervals.
                    continue;
                }
                for r in tree.find(iv.0..iv.1) {
                    hit = true;
                    let entry = r.data();
                    let gene_no = (*entry).0;
                    let nh = read.aux(b"NH");
                    let nh = nh.map_or(1, |aux| match aux {
                        Aux::I8(v) => v as u32,
                        Aux::U8(v) => v as u32,
                        Aux::I16(v) => v as u32,
                        Aux::U16(v) => v as u32,
                        Aux::I32(v) => v as u32,
                        Aux::U32(v) => v,
                        _ => 1,
                    });
                    if nh == 1 {
                        gene_nos_seen.insert(gene_no);
                    } else {
                        let hs = multimapper_dedup
                            .entry(gene_no)
                            .or_insert_with(HashSet::new);
                        hs.insert(read.qname().to_vec());
                    }
 
                    if each_read_counts_once {
                        break; // enable this (part 1 of 2) for each read hitting only once
                    }
                }
                if each_read_counts_once {
                    if hit {
                        //enable this (part 2 of 2) for each read hitting only once
                        break;
                    }
                }
            }
        }
        if !hit && !skipped {
            outside_count += 1;
        }
        for gene_no in gene_nos_seen.iter() {
            result[*gene_no as usize] += 1;
        }
    }
    for (gene_no, hs) in multimapper_dedup.iter() {
        result[*gene_no as usize] += hs.len() as u32;
    }
    Ok((result, outside_count))
}

///Injection point to do something with each read that actually matched
trait ReadCatcher {
    fn catch(&mut self, read: &bam::Record);
}

impl ReadCatcher for () {
    fn catch(&mut self, _read: &bam::Record) {}
}

///
/// count_reads_in_region_stranded
//
/// counts the stranded reads in a region,
/// matching them to the tree entries.
/// returns two vectors to be translated from gene_no
/// to gene_id: matching, reverse matching
fn count_reads_in_region_stranded<T>(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    gene_count: u32,
    each_read_counts_once: bool,
    mut read_catcher: T,
) -> Result<(Vec<u32>, Vec<u32>, u32)>
where
    T: ReadCatcher,
{
    let mut result_forward = vec![0; gene_count as usize];
    let mut result_reverse = vec![0; gene_count as usize];
    let mut multimapper_dedup_forward: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut multimapper_dedup_reverse: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut gene_nos_seen_forward = HashSet::<u32>::new();
    let mut gene_nos_seen_reverse = HashSet::<u32>::new();
    let mut outside_count = 0;
    let mut read: bam::Record = bam::Record::new();
    bam.fetch((tid, start as u64, stop as u64))?;
    while let Some(bam_result) = bam.read(&mut read) {
        bam_result?;
        // do not count multiple blocks matching in one gene multiple times
        gene_nos_seen_forward.clear();
        gene_nos_seen_reverse.clear();
        let mut hit = false;
        let mut skipped = false; //skipped are reads that don't belong to this region, but do get fetched by fetch.
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            skipped = true;
        }
        if !skipped {
            let mut gene_hit = false;
            let blocks = read.blocks();
            for iv in blocks.iter() {
                for r in tree.find(iv.0..iv.1) {
                    hit = true;
                    let entry = r.data();
                    let gene_no = (*entry).0;
                    let strand = (*entry).1; // this is 1 or -1
                    let nh = read.aux(b"NH");
                    let nh = nh.map_or(1, |aux| match aux {
                        Aux::I8(v) => v as u32,
                        Aux::U8(v) => v as u32,
                        Aux::I16(v) => v as u32,
                        Aux::U16(v) => v as u32,
                        Aux::I32(v) => v as u32,
                        Aux::U32(v) => v,
                        _ => 1,
                    });
                    if ((strand == 1) && !read.is_reverse()) || ((strand != 1) && read.is_reverse())
                    {
                        // read is in correct orientation - a hit
                        if nh == 1 {
                            gene_nos_seen_forward.insert(gene_no);
                            gene_hit = true;
                        } else {
                            let hs = multimapper_dedup_forward
                                .entry(gene_no)
                                .or_insert_with(HashSet::new);
                            hs.insert(read.qname().to_vec());
                            gene_hit = true;
                        }
                    } else {
                        // read is in incorrect orientation - not a hit
                        if nh == 1 {
                            gene_nos_seen_reverse.insert(gene_no);
                        } else {
                            let hs = multimapper_dedup_reverse
                                .entry(gene_no)
                                .or_insert_with(HashSet::new);
                            hs.insert(read.qname().to_vec());
                        }
                    }
                    if each_read_counts_once {
                        break; // enable this (part 1 of 2) for each read hitting only once
                    }
                }
                if each_read_counts_once {
                    if hit {
                        //enable this (part 2 of 2) for each read hitting only once
                        break;
                    }
                }
            }
            if gene_hit {
                read_catcher.catch(&read);
            }
        }

        if !hit && !skipped {
            outside_count += 1;
        }
        for gene_no in gene_nos_seen_forward.iter() {
            result_forward[*gene_no as usize] += 1;
        }
        for gene_no in gene_nos_seen_reverse.iter() {
            result_reverse[*gene_no as usize] += 1;
        }
    }
    for (gene_no, hs) in multimapper_dedup_forward.iter() {
        result_forward[*gene_no as usize] += hs.len() as u32;
    }
    for (gene_no, hs) in multimapper_dedup_reverse.iter() {
        result_reverse[*gene_no as usize] += hs.len() as u32;
    }
    Ok((result_forward, result_reverse, outside_count))
}


