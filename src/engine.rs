use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
mod chunked_genome;

use crate::bam_ext::BamRecordExtensions;
use crate::filters::ReadFilter;
use crate::gtf::Strand;
use crate::quantification::{Quant, Quantification};
use anyhow::{bail, Context, Result};
use bio::data_structures::interval_tree::IntervalTree;
use chunked_genome::{Chunk, ChunkedGenome};
use itertools::izip;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_htslib::bam::{self, Read};

use crate::gtf::GTFEntrys;
pub type OurTree = IntervalTree<u32, (u32, Strand)>;

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
pub struct IntervalResult {
    pub counts: HashMap<String, (f64, f64)>,
    pub filtered: usize,
    pub total: usize,
    pub outside: usize,
}

impl IntervalResult {
    fn new() -> Self {
        IntervalResult {
            counts: HashMap::new(),
            filtered: 0,
            total: 0,
            outside: 0,
        }
    }
}

pub struct Engine {
    reference_to_count_trees: HashMap<String, (OurTree, Vec<String>)>,
    reference_to_aggregation_trees: HashMap<String, (OurTree, Vec<String>)>,
    filters: Vec<crate::filters::Filter>,
    quantifier: Quantification,
}

impl Engine {
    pub fn from_gtf(
        mut gtf_entries: HashMap<String, GTFEntrys>,
        entry_kind: &str,
        entry_id_attribute: &str,
        aggregation_kind: &str,
        aggregation_id_attribute: &str,
        filters: Vec<crate::filters::Filter>,
        quantifier: Quantification,
    ) -> Result<Self> {
        let split_entries = gtf_entries
                .remove(aggregation_kind)
                .with_context(||format!("No GTF entries found for split feature {}. Necessary to know where to split the genome productivly.", aggregation_kind))?;
        let split_trees = build_trees_from_gtf(aggregation_id_attribute, split_entries)
            .context("Failed to build split trees from GTF")?;

        let feature_trees = if entry_kind == aggregation_kind {
            split_trees.clone()
        } else {
            build_trees_from_gtf(
                entry_id_attribute,
                gtf_entries.remove(entry_kind).expect("unreachable"),
            )
            .context("Failed to build feature trees")?
        };

        Ok(Engine {
            reference_to_count_trees: feature_trees,
            reference_to_aggregation_trees: split_trees,
            filters,
            quantifier,
        })
    }

    pub fn quantify_bam(
        &self,
        bam_path: impl AsRef<Path>,
        index_path: Option<&Path>,
        output_path: impl AsRef<Path>,
        write_output_bam: bool,
    ) -> Result<IntervalResult> {
        //check whether the bam file can be openend
        //and we need it for the chunking
        let bam_filename = bam_path.as_ref();
        let index_filename: Option<&Path> = index_path;
        let bam = crate::io::open_indexed_bam(bam_filename, index_filename.as_ref())?;
        let output_prefix: &Path = output_path.as_ref();
        let output_bam_info = match write_output_bam {
            true => {
                //create the output directory if it does not exist
                let od = output_prefix.join("annotated.bam.temp");
                if od.exists() {
                    //remove the whole tree
                    ex::fs::remove_dir_all(&od).with_context(|| {
                        format!("Failed to remove existing directory: {}", od.display())
                    })?;
                }
                std::fs::create_dir_all(&od).with_context(|| {
                    format!("Failed to create bam output directory: {}", od.display())
                })?;
                let header = rust_htslib::bam::Header::from_template(bam.header());

                Some((od, header))
            }

            false => None,
        };

        //perform the counting
        let cg = ChunkedGenome::new(&self.reference_to_aggregation_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
        let it: Vec<Chunk> = cg.iter().collect();
        let chunk_names = it.iter().map(|c| c.str_id()).collect::<Vec<_>>();
        let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
        let aggregated = pool.install(|| {
            let result: IntervalResult = it
                .into_par_iter()
                .map(|chunk| -> IntervalResult {
                    let mut bam =
                        crate::io::open_indexed_bam(bam_filename, index_filename).unwrap();
                    let (tree, gene_ids) = self.reference_to_count_trees.get(&chunk.chr).unwrap();
                    let counts = self.quantify_in_region(
                        &mut bam,
                        tree,
                        &gene_ids,
                        chunk.tid,
                        chunk.start,
                        chunk.stop,
                        output_bam_info
                            .as_ref()
                            .map(|(p, header)| (p.join(chunk.str_id().as_str()), header)),
                    );
                    counts.expect("Failure in quantification")
                })
                .reduce(IntervalResult::new, Self::aggregate);
            result
        });

        if let Some((temp_dir, output_header)) = output_bam_info {
            self.combine_temporary_bams(&chunk_names, temp_dir, output_prefix, output_header)?;
        }

        Ok(aggregated)
    }

    //called once per chunk.
    fn quantify_in_region(
        &self,
        bam: &mut rust_htslib::bam::IndexedReader,
        tree: &OurTree,
        gene_ids: &[String],
        tid: u32,
        start: u32,
        stop: u32,
        bam_output_path: Option<(PathBuf, &rust_htslib::bam::Header)>,
    ) -> Result<IntervalResult> {
        let mut result = vec![(0.0, 0.0); gene_ids.len()];
        let mut gene_nos_seen_match = HashSet::<u32>::new();
        let mut gene_nos_seen_reverse = HashSet::<u32>::new();
        let mut filtered_count = 0;
        let mut outside_count = 0;
        let mut total_count = 0;
        let mut read: bam::Record = bam::Record::new();
        let mut bam_out: Option<rust_htslib::bam::Writer> =
            bam_output_path.map(|(path, header)| {
                let mut writer = rust_htslib::bam::Writer::from_path(
                    path,
                    header,
                    rust_htslib::bam::Format::Bam,
                )
                .expect("Failed to create BAM writer");
                writer.set_threads(1).expect("Failed to set threads");
                writer
            });
        let mut quant = self.quantifier.clone();
        bam.fetch((tid, start as u64, stop as u64))?;
        'outer: while let Some(bam_result) = bam.read(&mut read) {
            bam_result?;
            total_count += 1;
            // do not count multiple blocks matching in one gene multiple times
            gene_nos_seen_match.clear();
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
                        let strand = entry.1;
                        let target = match (read.is_reverse(), strand) {
                            (_, Strand::Unstranded) => &mut gene_nos_seen_match,
                            (false, Strand::Plus) => &mut gene_nos_seen_match,
                            (false, Strand::Minus) => &mut gene_nos_seen_reverse,
                            (true, Strand::Plus) => &mut gene_nos_seen_reverse,
                            (true, Strand::Minus) => &mut gene_nos_seen_match,
                        };
                        target.insert(gene_no);
                    }
                }
            }
            if !hit && !skipped {
                outside_count += 1;
            }

            let weight = Self::decide_on_count(
                &read,
                &gene_nos_seen_match,
                &gene_nos_seen_reverse,
                &mut filtered_count,
                &self.filters,
                &mut quant,
            );
            if let Some(bam_out) = &mut bam_out {
                if !skipped {
                    let mut tag = String::new();
                    let mut first = true;
                    for (gene_no, w) in weight.0.iter() {
                        if !first {
                            tag.push(',');
                        } else {
                            first = false;
                        }
                        tag.push_str(&format!("{}={},", gene_ids[*gene_no as usize], w));
                    }
                    read.push_aux(b"XQ", rust_htslib::bam::record::Aux::String(&tag))
                        .expect("failed to push XQ tag");

                    let mut tag = String::new();
                    let mut first = true;
                    for (gene_no, w) in weight.1.iter() {
                        if !first {
                            tag.push(',');
                        } else {
                            first = false;
                        }
                        tag.push_str(&format!("{}={}", gene_ids[*gene_no as usize], w));
                    }
                    read.push_aux(b"XR", rust_htslib::bam::record::Aux::String(&tag))
                        .expect("failed to push XR tag");

                    bam_out.write(&read).expect("Failed to write BAM record");
                }
            } else if skipped {
                filtered_count += 1;
                continue 'outer; // skip this read
            }
            for (gene_id, w) in weight.0 {
                result[gene_id as usize].0 += w;
            }

            for (gene_id, w) in weight.1 {
                result[gene_id as usize].1 += w;
            }
        }

        Ok(IntervalResult {
            counts: gene_ids
                .iter()
                .enumerate()
                .map(|(idx, id)| (id.clone(), result[idx]))
                .collect(),
            filtered: filtered_count,
            total: total_count,
            outside: outside_count,
        })
    }

    fn decide_on_count<'a>(
        read: &bam::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
        filtered_count: &mut usize,
        filters: &Vec<crate::filters::Filter>,
        quantifier: &mut Quantification,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        let mut filtered = false;
        for f in filters.iter() {
            if f.remove_read(&read) {
                // if the read does not pass the filter, skip it
                *filtered_count += 1;
                filtered = true;
                break;
            }
        }
        if filtered {
            (
                gene_nos_seen_match.iter().map(|&id| (id, 0.0)).collect(),
                gene_nos_seen_reverse.iter().map(|&id| (id, 0.0)).collect(),
            )
        } else {
            quantifier.weight_read(read, gene_nos_seen_match, gene_nos_seen_reverse)
        }
    }

    fn aggregate(a: IntervalResult, b: IntervalResult) -> IntervalResult {
        IntervalResult {
            counts: add_dual_hashmaps(a.counts, b.counts),
            filtered: a.filtered + b.filtered,
            total: a.total + b.total,
            outside: a.outside + b.outside,
        }
    }

    fn combine_temporary_bams(
        &self,
        chunk_names: &[String],
        temp_dir: PathBuf,
        output_prefix: &Path,
        output_header: rust_htslib::bam::Header,
    ) -> Result<()> {
        //write the bam file from the generated chunks.
        // We need to do them in the right order. That means
        // we need to split the chunk_names by their reference,
        // iterate them in the order of the reference in the header
        // (and by position)
        let sorted_chunk_names = {
            let mut res = Vec::new();
            let hv: rust_htslib::bam::HeaderView =
                rust_htslib::bam::HeaderView::from_header(&output_header);
            for target in hv.target_names() {
                let target_name = String::from_utf8(target.to_vec())
                    .expect("Failed to convert target name to string");
                let target_chunks = chunk_names.iter().filter(|&x| {
                    x.split_once(":").expect("chunk name format unexpected").0 == target_name
                });
                //now verify that they're sorted by position
                let mut last_pos = 0;
                for chunk_name in target_chunks {
                    let parts: Vec<&str> = chunk_name.split(':').collect();
                    if parts.len() != 3 {
                        bail!(
                            "Chunk name {} is not in the expected format <chr>:<start>:<stop>",
                            chunk_name
                        );
                    }
                    let pos: u32 = parts[1]
                        .parse()
                        .context("Failed to parse position from chunk name")?;
                    if pos < last_pos {
                        bail!(
                            "Chunk names are not sorted by position: {} < {}",
                            pos,
                            last_pos
                        );
                    }
                    last_pos = pos;
                    res.push(chunk_name.clone())
                }
            }
            res
        };
        {
            let mut writer = rust_htslib::bam::Writer::from_path(
                output_prefix.join("annotated.bam"),
                &output_header,
                rust_htslib::bam::Format::Bam,
            )
            .expect("Failed to create BAM writer");
            for chunk_name in sorted_chunk_names {
                let mut input_bam =
                    rust_htslib::bam::Reader::from_path(temp_dir.join(chunk_name.as_str()))
                        .context("Failed to open a chunk of the output BAM to aggregate")?;
                let mut read = bam::Record::new();
                while let Some(bam_result) = input_bam.read(&mut read) {
                    bam_result.context("Failed to read BAM record")?;
                    writer.write(&read).context("Failed to write BAM record")?;
                }
            }
        }
        rust_htslib::bam::index::build(
            output_prefix.join("annotated.bam"),
            None,
            rust_htslib::bam::index::Type::Bai,
            0,
        )
        .context("Failed to build BAM index")?;
        //remove the temporary directory
        ex::fs::remove_dir_all(&temp_dir).with_context(|| {
            format!(
                "Failed to remove temporary BAM directory: {}",
                temp_dir.display()
            )
        })?;
        Ok(())
    }
}

fn add_dual_hashmaps(
    mut a: HashMap<String, (f64, f64)>,
    b: HashMap<String, (f64, f64)>,
) -> HashMap<String, (f64, f64)> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert((0., 0.));
        *x = (x.0 + v.0, x.1 + v.1);
    }
    a
}
