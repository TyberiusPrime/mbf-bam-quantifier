use super::{build_trees_from_gtf, IntervalCounter, IntervalResult, OurTree, Quant};
use crate::{
    bam_ext::BamRecordExtensions,
    config::{Input, Output},
    filters::{Filter, ReadFilter},
    quantification::IntervalIntermediateResult,
};
use anyhow::{bail, Context, Result};
use rust_htslib::bam::{self, Read};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    io::Write,
    path::Path,
};

fn default_split_feature() -> String {
    "gene".to_string()
}

fn default_split_id_attribute() -> String {
    "gene_id".to_string()
}

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Quantification {
    feature: String,
    id_attribute: String,
    #[serde(default = "default_split_feature")]
    split_feature: String,
    #[serde(default = "default_split_id_attribute")]
    split_id_attribute: String,
}

impl Quantification {
    fn open_output(&self, output: &Output) -> anyhow::Result<std::io::BufWriter<ex::fs::File>> {
        let output_file = format!("{}/counts.tsv", output.directory);
        ex::fs::create_dir_all(&output.directory)?;

        let output_file = ex::fs::File::create(&output_file)?;
        let mut out_buffer = std::io::BufWriter::new(output_file);
        out_buffer
            .write_all(format!("{}\tcount_correct\tcount_reverse\n", self.id_attribute).as_bytes())
            .expect("Failed to write header to output file");
        Ok(out_buffer)
    }
}

impl Quant for Quantification {
    fn quantify(
        &mut self,
        input: &Input,
        filters: &[Filter],
        output: &Output,
    ) -> anyhow::Result<()> {
        let mut gtf_entrys =
            input.read_gtf(&[self.feature.as_str(), self.split_feature.as_str()])?;

        let strand_info: Vec<_> = {
            let entries = gtf_entrys
                .get(&self.feature)
                .context("no GTF entries found for the specified feature")?;
            let ids = entries
                .vec_attributes
                .get(&self.id_attribute)
                .context("Missing id attribute")?
                .iter();
            entries
                .strand
                .iter()
                .zip(ids)
                .map(|(strand, id)| (id.to_string(), *strand))
                .collect()
        };

        let split_entries = gtf_entrys
                .remove(self.split_feature.as_str())
                .with_context(||format!("No GTF entries found for split feature {}. Necessary to know where to split the genome productivly.", &self.split_feature))?;
        let split_trees = build_trees_from_gtf(&self.split_id_attribute, split_entries)
            .context("Failed to build split trees from GTF")?;

        let feature_trees = if self.feature == self.split_feature {
            split_trees.clone()
        } else {
            build_trees_from_gtf(
                &self.id_attribute,
                gtf_entrys.remove(&self.feature).expect("unreachable"),
            )
            .context("Failed to build feature trees")?
        };
        let mut counts = StrandedCounter::count(
            Path::new(&input.bam),
            None,
            &feature_trees,
            split_trees,
            filters,
        )?; //these are counts in (forward, reverse) orientation.
            //We now need to to swap / sum them depending on what's in the gtf.

        for (entry_id, strand) in strand_info {
            match counts.counts.entry(entry_id) {
                std::collections::hash_map::Entry::Occupied(mut e) => match strand {
                    1 => {}
                    -1 => {
                        let c = e.get_mut();
                        *c = (c.1, c.0); // swap forward and reverse counts
                    }
                    0 => {
                        let c = e.get_mut();
                        *c = (c.0 + c.1, 0); // sum counts, set reverse to 0
                    }
                    _ => {
                        bail!("Unexpected strand value in GTF: {}", strand);
                    }
                },
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert((0, 0)); // initialize with zero counts
                }
            }
        }
        let mut out_buffer = self.open_output(output)?;

        let sorted_keys = {
            let mut keys: Vec<String> = counts.counts.keys().cloned().collect();
            keys.sort();
            keys
        };

        for feature_id in &sorted_keys {
            let count = counts.counts.get(feature_id).unwrap();
            out_buffer
                .write_all(format!("{}\t{}\t{}\n", feature_id, count.0, count.1).as_bytes())?;
        }
        //todo: write total / outside to some sidechannel?

        Ok(())
    }

    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.input.gtf.is_none() {
            bail!("GTF file is required for featurecount quantification");
        }

        Ok(())
    }
}

/// the easiest counter possible.
/// For every bam entry, figure out which intervals it covers,
/// then add 1 to the count of interval
/// either to the .0 (forward oriented read) or .1 (reverse oriented read)
/// We add them all up, and afterwards swap / aggregate depending on
/// what the gtf tells us...
#[derive(Debug)]
struct StrandedCounter {}

impl IntervalCounter for StrandedCounter {
    type OutputType = (u32, u32);

    #[allow(clippy::nonminimal_bool)]
    fn count_in_region(
        bam: &mut rust_htslib::bam::IndexedReader,
        tree: &OurTree,
        tid: u32,
        start: u32,
        stop: u32,
        gene_ids_len: u32, //how many are in the tree
        filters: &[crate::filters::Filter],
    ) -> Result<IntervalIntermediateResult<Self::OutputType>> {
        let mut result = vec![(0, 0); gene_ids_len as usize];
        let mut gene_nos_seen_forward = HashSet::<u32>::new();
        let mut gene_nos_seen_reverse = HashSet::<u32>::new();
        let mut outside_count = 0;
        let mut total_count = 0;
        let mut read: bam::Record = bam::Record::new();
        bam.fetch((tid, start as u64, stop as u64))?;
        'outer: while let Some(bam_result) = bam.read(&mut read) {
            bam_result?;
            for f in filters.iter() {
                if f.remove_read(&read) {
                    // if the read does not pass the filter, skip it
                    continue 'outer;
                }
            }
            // do not count multiple blocks matching in one gene multiple times
            gene_nos_seen_forward.clear();
            gene_nos_seen_reverse.clear();

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
                        let gene_no = entry.0;
                        if read.is_reverse() {
                            gene_nos_seen_reverse.insert(gene_no);
                        } else {
                            gene_nos_seen_forward.insert(gene_no);
                        }
                    }
                }
            }
            if !hit && !skipped {
                outside_count += 1;
            }
            total_count += 1;
            for gene_no in gene_nos_seen_forward.iter() {
                result[*gene_no as usize].0 += 1;
            }
            for gene_no in gene_nos_seen_reverse.iter() {
                result[*gene_no as usize].1 += 1;
            }
        }
        Ok(IntervalIntermediateResult {
            counts: result,
            total: total_count,
            outside: outside_count,
        })
    }

    fn aggregate(
        a: IntervalResult<Self::OutputType>,
        b: IntervalResult<Self::OutputType>,
    ) -> IntervalResult<Self::OutputType> {
        IntervalResult {
            counts: add_dual_hashmaps(a.counts, b.counts),
            total: a.total + b.total,
            outside: a.outside + b.outside,
        }
    }
}

fn add_dual_hashmaps(
    mut a: HashMap<String, (u32, u32)>,
    b: HashMap<String, (u32, u32)>,
) -> HashMap<String, (u32, u32)> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert((0, 0));
        *x = (x.0 + v.0, x.1 + v.1);
    }
    a
}
