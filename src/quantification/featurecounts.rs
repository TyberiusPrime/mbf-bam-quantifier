mod chunked_genome;
mod counters;

use itertools::izip;
use std::collections::HashMap;

use super::Quant;
use crate::config::{Input, Output};
use crate::gtf::GTFEntrys;
use anyhow::{bail, Context, Result};
use bio::data_structures::interval_tree::IntervalTree;
use serde::{Deserialize, Serialize};
use std::io::Write;

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Quantification {
    feature: String,
    id_attribute: String,
}

impl Quantification {
    fn open_output(&self, output: &Output) -> anyhow::Result<std::io::BufWriter<ex::fs::File>> {
        let output_file = format!("{}/counts.tsv", output.directory);
        ex::fs::create_dir_all(&output.directory)?;

        let output_file = ex::fs::File::create(&output_file)?;
        let mut out_buffer = std::io::BufWriter::new(output_file);
        out_buffer
            .write_all(b"reference\tcount\n")
            .expect("Failed to write header to output file");
        Ok(out_buffer)
    }
}

impl Quant for Quantification {
    fn quantify(&mut self, input: &Input, output: &Output) -> anyhow::Result<()> {
        // Implement the quantification logic here
        let bam_with_index = input.get_indexed_bam_reader()?;

        let mut gtf_entrys = if self.feature == "gene" {
            input.read_gtf(&[self.feature.as_str()])?
        } else {
            input.read_gtf(&["gene", &self.feature])?
        };

        let gene_trees =
            build_trees_from_gtf(&self.id_attribute, gtf_entrys.remove("gene").context("No GTF entries found for `gene`s. Necessary to know where to split the genome productivly.")?)
                .context("Failed to build trees from GTF")?;
        let feature_trees = if self.feature == "gene" {
            gene_trees.clone()
        } else {
            build_trees_from_gtf(
                &self.id_attribute,
                gtf_entrys
                    .remove(&self.feature)
                    .context("No GTF entries found for the specified feature")?,
            )?
        };

        let counts = counters::count_reads_stranded(
            &input.bam,
            None,
            gene_trees,
            feature_trees,
            true,
            None,
        );


        Ok(())
    }

    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.input.gtf.is_none() {
            bail!("GTF file is required for featurecount quantification");
        }

        Ok(())
    }
}

// OurTree stores an interval tree
/// a gene_no (ie. an index into a vector of gene_ids)
/// and the strand (+1/ -1, 0)
pub type OurTree = IntervalTree<u32, (u32, i8)>;

pub fn build_trees_from_gtf(
    id_attribute: &str,
    gtf_entries: GTFEntrys,
) -> Result<HashMap<String, (OurTree, Vec<String>)>> {
    let mut trees: HashMap<u32, OurTree> = HashMap::new();
    let mut gene_nos_by_chr = HashMap::new();
    for (line, (seq_name_cat_id, gene_id, start, end, strand)) in izip!(
        gtf_entries.seqname.values.iter(),
        gtf_entries
            .vec_attributes
            .get(id_attribute)
            .context("Missing id attribute")?
            .iter(),
        gtf_entries.start.iter(),
        gtf_entries.end.iter(),
        gtf_entries.strand.iter(),
    )
    .enumerate()
    {
        let (gene_nos, genes_in_order) = gene_nos_by_chr
            .entry(*seq_name_cat_id)
            .or_insert_with(|| (HashMap::new(), Vec::new()));
        let gene_no = gene_nos.entry(gene_id).or_insert_with(|| {
            genes_in_order.push(gene_id.clone());
            genes_in_order.len() - 1
        });

        let tree = trees
            .entry(*seq_name_cat_id)
            .or_insert_with(|| OurTree::new());
        let start: u32 = (*start)
            .try_into()
            .context("Start value is not a valid u64")?;
        let end: u32 = (*end).try_into().context("End value is not a valid u64")?;

        tree.insert(
            start..end, //these are already 0-based
            (*gene_no as u32, *strand),
        )
    }
    let res: HashMap<String, _> = trees
        .into_iter()
        .map(|(k, v)| (gtf_entries.seqname.cat_from_value(k), v))
        .zip(gene_nos_by_chr.into_iter().map(|(_k, (gene_nos, genes_in_order))| {
            genes_in_order
        }))
        .map(|((seq_name, tree), genes_in_order)| {
            (seq_name, (tree, genes_in_order))
        })
        .collect();

    Ok((res))
}

fn add_hashmaps(mut a: HashMap<String, u32>, b: HashMap<String, u32>) -> HashMap<String, u32> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert(0);
        *x += v;
    }
    a
}
