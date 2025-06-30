//mod featurecounts;
mod chunked_genome;
mod references;
mod stranded;
mod unstranded;

use crate::config::{Config, Input, Output};
use crate::gtf::GTFEntrys;
use anyhow::{Context, Result};
use bio::data_structures::interval_tree::IntervalTree;
use chunked_genome::{Chunk, ChunkedGenome};
use enum_dispatch::enum_dispatch;
use itertools::izip;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use std::path::Path;

#[enum_dispatch(Quantification)]
pub trait Quant {
    fn quantify(
        &mut self,
        input: &Input,
        filters: &[crate::filters::Filter],
        output: &Output,
    ) -> anyhow::Result<()>;
    fn check(&self, _config: &Config) -> anyhow::Result<()> {
        Ok(())
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, serde::Serialize)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Quantification {
    #[serde(alias = "references")]
    References(references::Quantification),
    /* #[serde(alias = "featurecounts")]
    FeatureCounts(featurecounts::Quantification), */
    #[serde(alias = "gtf_unstranded")]
    GTFUnstranded(unstranded::Quantification),

    #[serde(alias = "gtf_stranded")]
    GTFStranded(stranded::Quantification),
}

impl Quantification {}

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
    for (seq_name_cat_id, gene_id, start, end, strand) in izip!(
        gtf_entries.seqname.values.iter(),
        gtf_entries
            .vec_attributes
            .get(id_attribute)
            .context("Missing id attribute")?
            .iter(),
        gtf_entries.start.iter(),
        gtf_entries.end.iter(),
        gtf_entries.strand.iter(),
    ) {
        println!(
            "{} {} {}",
            *seq_name_cat_id,
            gene_id,
            gtf_entries.seqname.cat_from_value(*seq_name_cat_id)
        );
        let (gene_nos, genes_in_order) = gene_nos_by_chr
            .entry(*seq_name_cat_id)
            .or_insert_with(|| (HashMap::new(), Vec::new()));
        let gene_no = gene_nos.entry(gene_id).or_insert_with(|| {
            genes_in_order.push(gene_id.clone());
            genes_in_order.len() - 1
        });

        let tree = trees.entry(*seq_name_cat_id).or_default();
        let start: u32 = (*start)
            .try_into()
            .context("Start value is not a valid u64")?;
        let end: u32 = (*end).try_into().context("End value is not a valid u64")?;

        tree.insert(
            start..end, //these are already 0-based
            (*gene_no as u32, *strand),
        )
    }
    let res: Result<HashMap<String, _>> = trees
        .into_iter()
        .map(|(seq_name_cat_id, tree)| {
            let seq_name = gtf_entries.seqname.cat_from_value(seq_name_cat_id);
            let gene_nos = gene_nos_by_chr
                .remove(&seq_name_cat_id)
                .context("Missing gene numbers for sequence name")?;
            Ok((seq_name, (tree, gene_nos.1)))
        })
        .collect();

    res
}

#[derive(Debug)]
struct IntervalResult<T: std::marker::Send + std::marker::Sync + std::fmt::Debug> {
    counts: HashMap<String, T>,
    total: usize,
    outside: usize,
}

impl<T: std::marker::Send + std::marker::Sync + std::fmt::Debug> IntervalResult<T> {
    fn new() -> Self {
        IntervalResult {
            counts: HashMap::new(),
            total: 0,
            outside: 0,
        }
    }
}

#[derive(Debug)]
struct IntervalIntermediateResult<T: std::marker::Send + std::marker::Sync + std::fmt::Debug> {
    counts: Vec<T>,
    total: usize,
    outside: usize,
}

trait IntervalCounter
where
    Self::OutputType: std::marker::Send + std::marker::Sync + std::fmt::Debug,
{
    type OutputType;

    fn count(
        bam_filename: &Path,
        index_filename: Option<&Path>,
        interval_trees: &HashMap<String, (OurTree, Vec<String>)>,
        split_trees: HashMap<String, (OurTree, Vec<String>)>,
        filters: &[crate::filters::Filter],
    ) -> Result<IntervalResult<Self::OutputType>> {
        //check whether the bam file can be openend
        //and we need it for the chunking
        let bam = crate::io::open_indexed_bam(bam_filename, index_filename)?;

        //perform the counting
        let cg = ChunkedGenome::new(split_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
        let it: Vec<Chunk> = cg.iter().collect();
        let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
        let aggregated = pool.install(|| {
            let result: IntervalResult<Self::OutputType> = it
                .into_par_iter()
                .map(|chunk| -> IntervalResult<Self::OutputType> {
                    let mut bam =
                        crate::io::open_indexed_bam(bam_filename, index_filename).unwrap();
                    let (tree, gene_ids) = interval_trees.get(&chunk.chr).unwrap();

                    let local_result: IntervalIntermediateResult<Self::OutputType> =
                        Self::count_in_region(
                            &mut bam,
                            //&chunk.tree,
                            tree,
                            chunk.tid,
                            chunk.start,
                            chunk.stop,
                            gene_ids.len() as u32,
                            filters,
                        )
                        .expect("Failure during read counting");
                    let result = local_result
                        .counts
                        .into_iter()
                        .enumerate()
                        .map(|(k, v)| (gene_ids[k].clone(), v))
                        .collect::<HashMap<String, Self::OutputType>>();
                    IntervalResult {
                        counts: result,
                        total: local_result.total,
                        outside: local_result.outside,
                    }
                })
                //todo: can we convert this into an error collector?
                .reduce(IntervalResult::new, Self::aggregate);
            Ok(result)
        });
        aggregated
    }

    fn count_in_region(
        bam: &mut rust_htslib::bam::IndexedReader,
        tree: &OurTree,
        tid: u32,
        start: u32,
        stop: u32,
        gene_ids_len: u32, //how many are in the tree
        filters: &[crate::filters::Filter],
    ) -> Result<IntervalIntermediateResult<Self::OutputType>>;

    fn aggregate(
        a: IntervalResult<Self::OutputType>,
        b: IntervalResult<Self::OutputType>,
    ) -> IntervalResult<Self::OutputType>;
}
