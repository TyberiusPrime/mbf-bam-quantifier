use super::OurTree;
use anyhow::Result;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::convert::TryFrom;

pub struct ChunkedGenome<'a> {
    trees: Option<&'a HashMap<String, (OurTree, Vec<String>)>>,
    bam: bam::IndexedReader,
    chromosomes: Vec<String>,
}

pub fn tids_with_reads(bam: &mut bam::IndexedReader) -> Result<Vec<u32>> {
    let keep_tids: Vec<u32> = bam
        .index_stats()?
        .iter()
        .filter_map(|(tid, _length, mapped_count, _unmapped_count)| {
            if *mapped_count > 0 {
                //since we're not dealing with the last one, which is -1 in i64
                Some((*tid).try_into().expect("SAM tid should fit into u32"))
            } else {
                None
            }
        })
        .collect();
    Ok(keep_tids)
}

impl<'a> ChunkedGenome<'a> {
    ///create a new chunked genome for iteration
    ///if you pass in a tree, it is guaranteed that the splits happen
    ///between entries of the tree, not inside.
    pub fn new(
        trees: &'a HashMap<String, (OurTree, Vec<String>)>,
        mut bam: bam::IndexedReader,
    ) -> Result<ChunkedGenome<'a>> {
        let keep_tids = tids_with_reads(&mut bam)?;
        let chrs_in_tree_and_bam = trees
            .keys()
            .filter(|x| {
                bam.header()
                    .tid(x.as_bytes())
                    .map(|tid| keep_tids.iter().any(|x| x == &tid))
                    .unwrap_or(false)
            })
            .cloned()
            .collect();
        Ok(ChunkedGenome {
            chromosomes: chrs_in_tree_and_bam,
            trees: Some(trees),
            bam,
        })
    }

    /* pub fn new_without_tree(bam: bam::IndexedReader) -> ChunkedGenome {
        ChunkedGenome {
            trees: None,
            chromosomes: bam
                .header()
                .target_names()
                .iter()
                .map(|x| str::from_utf8(x).unwrap().to_string())
                .collect(),
            bam,
        }
    } */

    pub fn iter(&self) -> ChunkedGenomeIterator {
        ChunkedGenomeIterator {
            cg: self,
            it: self.chromosomes.iter(),
            last_start: 0,
            last_tid: 0,
            last_chr_length: 0,
            last_chr: "".to_string(),
        }
    }
}

pub struct ChunkedGenomeIterator<'a> {
    cg: &'a ChunkedGenome<'a>,
    it: std::slice::Iter<'a, String>,
    last_start: u32,
    last_chr: String,
    last_tid: u32,
    last_chr_length: u32,
}
#[derive(Debug)]
pub struct Chunk {
    pub chr: String,
    pub tid: u32,
    pub start: u32,
    pub stop: u32,
}

impl Chunk {
    pub fn new(chr: String, tid: u32, start: u32, stop: u32) -> Chunk {
        Chunk {
            chr,
            tid,
            start,
            stop,
        }
    }
    pub fn str_id(&self) -> String {
        format!("{}:{}:{}", self.chr, self.start, self.stop)
    }

    pub fn interval_outside(&self, start: u32, stop:u32) -> bool {
        self.start > stop || self.stop < start
                /* if (iv.1 < chunk.start)
                    || iv.0 >= chunk.stop
                    || ((iv.0 < chunk.start) && (iv.1 >= chunk.start)) */
    }
}

impl Iterator for ChunkedGenomeIterator<'_> {
    type Item = Chunk;
    fn next(&mut self) -> Option<Chunk> {
        let chunk_size = 10_000_000;
        if self.last_start >= self.last_chr_length {
            let next_chr = self.it.next()?;
            let tid = self.cg.bam.header().tid(next_chr.as_bytes()).unwrap();
            let chr_length = u32::try_from(self.cg.bam.header().target_len(tid).unwrap())
                .expect("Not u64 contig length aware");
            self.last_tid = tid;
            self.last_chr_length = chr_length;
            self.last_chr = next_chr.to_string();
            self.last_start = 0;
        }

        let mut stop = self.last_start + chunk_size;
        if self.cg.trees.is_some() {
            let (next_tree, _next_gene_ids) =
                self.cg.trees.as_ref().unwrap().get(&self.last_chr).unwrap();
            loop {
                ////this has been adjusted not to cut genes in half
                //cut gene in half?
                //option 0 for that is to pass in the gene intervals as well
                //just for constructing the chunks
                let overlapping = next_tree.find(stop..stop + 1).next();
                match overlapping {
                    None => break,
                    Some(entry) => {
                        let iv = entry.interval();
                        if iv.end + 1 < stop {
                            panic!("WHAT?");
                        }
                        stop = iv.end + 1;
                    }
                }
            }
        }
        let c = Chunk {
            chr: self.last_chr.clone(),
            tid: self.last_tid,
            start: self.last_start,
            stop,
        };
        self.last_start = stop;
        Some(c)
    }
}
