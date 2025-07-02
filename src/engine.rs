use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
mod chunked_genome;

use crate::bam_ext::BamRecordExtensions;
use crate::extractors::UMIExtractor;
use crate::filters::ReadFilter;
use crate::gtf::Strand;
use crate::quantification::{Quant, Quantification};
use anyhow::{bail, Context, Result};
use bio::data_structures::interval_tree::IntervalTree;
use chunked_genome::{Chunk, ChunkedGenome};
use itertools::{izip, Itertools};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_htslib::bam::{self, Read};

use crate::gtf::GTFEntrys;
pub type OurTree = IntervalTree<u32, (u32, Strand)>;

pub fn build_trees_from_gtf(
    id_attribute: &str,
    gtf_entries: &GTFEntrys,
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

//first merge intervals by id_attribute, then build the tree
pub fn build_trees_from_gtf_merged(
    id_attribute: &str,
    gtf_entries: &GTFEntrys,
) -> Result<HashMap<String, (OurTree, Vec<String>)>> {
    let mut intervals: HashMap<u32, HashMap<String, (u32, u32)>> = HashMap::new();

    for (seq_name_cat_id, gene_id, start, end) in izip!(
        gtf_entries.seqname.values.iter(),
        gtf_entries
            .vec_attributes
            .get(id_attribute)
            .context("Missing id attribute")?
            .iter(),
        gtf_entries.start.iter(),
        gtf_entries.end.iter(),
    ) {
        let by_chr_interval = intervals
            .entry(*seq_name_cat_id)
            .or_insert_with(HashMap::new);
        let start: u32 = (*start)
            .try_into()
            .context("Start value is not a valid u32")?;
        let end: u32 = (*end)
            .try_into()
            .context("Start value is not a valid u32")?;

        match by_chr_interval.entry(gene_id.to_string()) {
            std::collections::hash_map::Entry::Occupied(mut e) => {
                let last = *e.get();
                let new = ((start).min(last.0), (end).max(last.1));
                *e.get_mut() = new;
            }
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert((start, end));
            }
        }
    }

    let mut res = HashMap::new();
    for (seq_name_cat_id, intervals) in intervals.iter() {
        let seq_name = gtf_entries.seqname.cat_from_value(*seq_name_cat_id);
        let mut tree = OurTree::new();
        for (start, stop) in intervals.values() {
            tree.insert(*start..*stop, (0, Strand::Unstranded));
        }
        res.insert(seq_name, (tree, vec!["ignored".to_string()]));
    }
    Ok(res)
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

pub enum AnnotatedRead {
    Filtered,
    NotInRegion,
    Counted(AnnotatedReadInfo),
    Duplicate,
}

pub struct AnnotatedReadInfo {
    pub corrected_position: u32, // clipping corrected position
    pub genes_hit_correct: HashMap<String, f32>, //todo: optimize!
    pub genes_hit_reverse: HashMap<String, f32>, //todo: optimize!
    pub umi: Option<String>,     // Optional: What's it's UMI.
    pub barcode: Option<String>, // Optional: What's it's cellbarcode
    pub mapping_priority: (u8, u8),
}

pub enum ReadToGeneMatcher {
    TreeMatcher(TreeMatcher),
    TagMatcher(TagMatcher),
}

impl ReadToGeneMatcher {
    fn generate_chunks(&self, bam: rust_htslib::bam::IndexedReader) -> Vec<Chunk> {
        match self {
            ReadToGeneMatcher::TreeMatcher(matcher) => matcher.generate_chunks(bam),
            ReadToGeneMatcher::TagMatcher(matcher) => matcher.generate_chunks(bam),
        }
    }

    fn hits(
        &self,
        chr: &str,
        chunk_start: u32,
        chunk_stop: u32,
        read: &rust_htslib::bam::record::Record,
        quant_is_reverse: bool,
    ) -> Result<(HashSet<String>, HashSet<String>)> {
        match self {
            ReadToGeneMatcher::TreeMatcher(matcher) => {
                matcher.hits(chr, chunk_start, chunk_stop, read, quant_is_reverse)
            }
            ReadToGeneMatcher::TagMatcher(matcher) => {
                matcher.hits(chr, chunk_start, chunk_stop, read, quant_is_reverse)
            }
        }
    }
}

pub struct Engine {
    quantifier: Quantification,
    filters: Vec<crate::filters::Filter>,
    matcher: ReadToGeneMatcher,
    umi_extractor: crate::extractors::UMIExtraction,
}

impl Engine {
    pub fn from_gtf(
        mut gtf_entries: HashMap<String, GTFEntrys>,
        entry_kind: &str,
        entry_id_attribute: &str,
        aggregation_id_attribute: &str,
        filters: Vec<crate::filters::Filter>,
        quantifier: Quantification,
        umi_extractor: crate::extractors::UMIExtraction,
    ) -> Result<Self> {
        let feature_entries =  gtf_entries
                .remove(entry_kind)
                .with_context(||format!("No GTF entries found for feature {}. Necessary to know where to split the genome productivly.", entry_kind))?;

        let feature_trees = build_trees_from_gtf(aggregation_id_attribute, &feature_entries)
            .context("Failed to build feature trees from GTF")?;

        let split_trees = if aggregation_id_attribute == entry_id_attribute {
            feature_trees.clone()
        } else {
            build_trees_from_gtf_merged(aggregation_id_attribute, &feature_entries)
                .context("Failed to build split trees from GTF")?
        };

        let feature_trees = if entry_id_attribute == aggregation_id_attribute {
            split_trees.clone()
        } else {
            build_trees_from_gtf(
                entry_id_attribute,
                &gtf_entries.remove(entry_kind).expect("unreachable"),
            )
            .context("Failed to build feature trees")?
        };

        Ok(Engine {
            matcher: ReadToGeneMatcher::TreeMatcher(TreeMatcher {
                reference_to_count_trees: feature_trees,
                reference_to_aggregation_trees: split_trees,
            }),
            filters,
            quantifier,
            umi_extractor,
        })
    }
    pub fn from_references(
        references: Vec<(String, u64)>,
        filters: Vec<crate::filters::Filter>,
        quantifier: Quantification,
        umi_extractor: crate::extractors::UMIExtraction,
    ) -> Result<Self> {
        let mut trees: HashMap<String, (OurTree, Vec<String>)> = HashMap::new();
        for (seq_name, length) in references.iter() {
            let mut tree = OurTree::new();
            let start = 0;
            let stop = *length as u32;
            tree.insert(
                start..stop,             //these are already 0-based
                (0, Strand::Unstranded), // no gene numbers here, just a dummy
            );
            let genes_in_order = vec![seq_name.clone()];
            trees.insert(seq_name.clone(), (tree, genes_in_order));
        }

        Ok(Engine {
            matcher: ReadToGeneMatcher::TreeMatcher(TreeMatcher {
                reference_to_count_trees: trees.clone(),
                reference_to_aggregation_trees: trees,
            }),
            filters,
            quantifier,
            umi_extractor,
        })
    }

    pub fn from_bam_tag(
        tag: [u8; 2],
        filters: Vec<crate::filters::Filter>,
        quantifier: Quantification,
        umi_extractor: crate::extractors::UMIExtraction,
    ) -> Self {
        assert!(
            tag[0] == tag[0].to_ascii_uppercase(),
            "BAM tag must be uppercase"
        );
        assert!(
            tag[1] == tag[1].to_ascii_uppercase(),
            "BAM tag must be uppercase"
        );
        Engine {
            matcher: ReadToGeneMatcher::TagMatcher(TagMatcher { tag }),
            filters,
            quantifier,
            umi_extractor,
        }
    }

    pub fn quantify_bam(
        &self,
        bam_path: impl AsRef<Path>,
        index_path: Option<&Path>,
        output_path: impl AsRef<Path>,
        write_output_bam: bool,
        max_skip_len: u32,
    ) -> Result<IntervalResult> {
        //check whether the bam file can be openend
        //and we need it for the chunking
        let bam_filename = bam_path.as_ref();
        let index_filename: Option<&Path> = index_path;
        let bam = crate::io::open_indexed_bam(bam_filename, index_filename.as_ref())?;
        // prepare output bam (if requested)
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
        let chunks = self.matcher.generate_chunks(bam);
        //perform the counting
        let chunk_names = chunks.iter().map(|c| c.str_id()).collect::<Vec<_>>();

        let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
        let aggregated = pool.install(|| {
            let result: IntervalResult = chunks
                .into_par_iter()
                .map(|chunk| -> IntervalResult {
                    let mut bam =
                        crate::io::open_indexed_bam(bam_filename, index_filename).unwrap();
                    let mut annotated_reads = self
                        .annotate_reads(
                            &mut bam,
                            &chunk.chr,
                            chunk.tid,
                            chunk.start,
                            chunk.stop,
                            max_skip_len,
                        )
                        .expect("Failure in quantification");

                    let filtered = annotated_reads
                        .iter()
                        .filter(|(read, _org_index)| matches!(read, AnnotatedRead::Filtered))
                        .count();
                    let total = annotated_reads.len();
                    let outside = annotated_reads
                        .iter()
                        .filter(|(read, _org_index)| matches!(read, AnnotatedRead::NotInRegion))
                        .count();

                    //annotated_reads.retain(|read| matches!(read, AnnotatedRead::Counted(_)));
                    annotated_reads.sort_by_key(|(read, org_index)| match read {
                        AnnotatedRead::Counted(info) => info.corrected_position as u64,
                        _ => u32::MAX as u64 + 1,
                    });

                    let mut last_pos = None;
                    let mut change_indices = Vec::new();
                    let mut found_filtered = false;
                    for (ii, (read, org_index)) in annotated_reads.iter().enumerate() {
                        let info = match read {
                            AnnotatedRead::Counted(info) => info,
                            _ => {
                                //we have reached the counted sector.
                                change_indices.push(ii);
                                found_filtered = true;
                                break;
                            }
                        };
                        let new_pos = Some(info.corrected_position);
                        if new_pos != last_pos {
                            change_indices.push(ii);
                            last_pos = Some(info.corrected_position)
                        }
                    }
                    if !found_filtered {
                        change_indices.push(annotated_reads.len());
                    }

                    for (start, stop) in change_indices.iter().tuple_windows() {
                        self.quantifier
                            .weight_read_group(&mut annotated_reads[*start..*stop])
                            .expect("weighting failed");
                    }

                    let mut counts: HashMap<String, (f64, f64)> = HashMap::new();
                    for (read, _org_index) in &annotated_reads {
                        match read {
                            AnnotatedRead::Duplicate
                            | AnnotatedRead::Filtered
                            | AnnotatedRead::NotInRegion => {}
                            AnnotatedRead::Counted(info) => {
                                for (gene, weight) in &info.genes_hit_correct {
                                    let entry =
                                        counts.entry(gene.to_string()).or_insert((0.0, 0.0));
                                    entry.0 += *weight as f64; // count forward
                                }
                                for (gene, weight) in &info.genes_hit_reverse {
                                    let entry =
                                        counts.entry(gene.to_string()).or_insert((0.0, 0.0));
                                    entry.1 += *weight as f64; // count reverse
                                }
                            }
                        }
                    }
                    if let Some((out_bam_path, header)) = output_bam_info.as_ref() {
                        Self::write_annotated_reads(
                            &mut bam,
                            &chunk,
                            &annotated_reads,
                            &out_bam_path,
                            header,
                            max_skip_len,
                        )
                        .expect("Failed to write output bam");
                    }
                    //todo: bam writing...
                    //
                    IntervalResult {
                        counts,
                        filtered,
                        total,
                        outside,
                    }
                })
                .reduce(IntervalResult::new, Self::aggregate);
            result
        });

        if let Some((temp_dir, output_header)) = output_bam_info {
            combine_temporary_bams(&chunk_names, temp_dir, output_prefix, output_header)?;
        }

        Ok(aggregated)
    }

    fn aggregate(a: IntervalResult, b: IntervalResult) -> IntervalResult {
        IntervalResult {
            counts: add_dual_hashmaps(a.counts, b.counts),
            filtered: a.filtered + b.filtered,
            total: a.total + b.total,
            outside: a.outside + b.outside,
        }
    }

    fn annotate_reads(
        &self,
        bam: &mut rust_htslib::bam::IndexedReader,
        chr: &str,
        tid: u32,
        start: u32,
        stop: u32,
        max_skip_len: u32,
    ) -> Result<Vec<(AnnotatedRead, usize)>> {
        let mut res = Vec::new();

        let mut read: bam::Record = bam::Record::new();
        let mut ii = 0;
        bam.fetch((
            tid,
            start as u64,
            (stop + max_skip_len) as u64, // this is a correctness issue
                                          // We'll read unnecessary reads here,
                                          // just to discard them because their correceted position is beyond the region
                                          // but there might be a read there that is.
        ))?;
        'outer: while let Some(bam_result) = bam.read(&mut read) {
            bam_result?;

            let corrected_position = read.corrected_pos(max_skip_len);

            if corrected_position >= 0 {
                if ((corrected_position as u32) < start) || ((corrected_position as u32) >= stop) {
                    //this ensures we count a read only in the chunk where it's left most
                    //pos is in.
                    res.push((AnnotatedRead::NotInRegion, ii));
                    continue;
                }

                for f in self.filters.iter() {
                    if f.remove_read(&read) {
                        // if the read does not pass the filter, skip it
                        res.push((AnnotatedRead::Filtered, ii));
                        continue 'outer;
                    }
                }

                let (genes_hit_correct, genes_hit_reverse) =
                    self.matcher
                        .hits(chr, start, stop, &read, self.quantifier.reverse())?;

                let info = AnnotatedReadInfo {
                    corrected_position: corrected_position as u32,
                    genes_hit_correct: genes_hit_correct.into_iter().map(|id| (id, 1.0)).collect(),
                    genes_hit_reverse: genes_hit_reverse.into_iter().map(|id| (id, 1.0)).collect(),
                    umi: self.umi_extractor.extract(&read),
                    barcode: None,
                    mapping_priority: (
                        read.no_of_mapping_coordinates().try_into().unwrap_or(255),
                        read.mapq(),
                    ),
                };
                res.push((AnnotatedRead::Counted(info), ii));
                ii += 1;
            } else {
                panic!("Did not expect to see an unaligned read in a chunk");
            }
        }
        Ok(res)
    }

    fn write_annotated_reads(
        bam: &mut rust_htslib::bam::IndexedReader,
        chunk: &Chunk,
        annotated_reads: &[(AnnotatedRead, usize)],
        out_bam_path: &Path,
        header: &rust_htslib::bam::Header,
        max_skip_len: u32,
    ) -> Result<()> {
        let mut out_bam = rust_htslib::bam::Writer::from_path(
            out_bam_path.join(chunk.str_id()),
            &header,
            rust_htslib::bam::Format::Bam,
        )?;
        let idx_to_annotated: HashMap<usize, &AnnotatedRead> = annotated_reads
            .iter()
            .map(|(read, org_index)| (*org_index, read))
            .collect();
        bam.fetch((
            chunk.tid,
            chunk.start as u64,
            (chunk.stop + max_skip_len) as u64, // this is a correctness issue
                                                // We'll read unnecessary reads here,
                                                // just to discard them because their correceted position is beyond the region
                                                // but there might be a read there that is.
        ))?;
        let mut read = bam::Record::new();
        let mut ii = 0;
        while let Some(bam_result) = bam.read(&mut read) {
            bam_result?;
            if let Some(&anno_read) = idx_to_annotated.get(&ii) {
                match anno_read {
                    AnnotatedRead::NotInRegion => continue,
                    AnnotatedRead::Filtered => {}
                    AnnotatedRead::Duplicate => {
                        read.push_aux(b"XD", rust_htslib::bam::record::Aux::U8(1))?;
                    }
                    AnnotatedRead::Counted(info) => {
                        //we have a read that was annotated
                        //write it to the output bam
                        let mut tag = String::new();
                        let mut first = true;
                        for (gene, weight) in &info.genes_hit_correct {
                            if !first {
                                tag.push(',')
                            }
                            first = false;
                            tag.push_str(&format!("{}={:.2}", gene, weight));
                        }
                        read.push_aux(b"XQ", rust_htslib::bam::record::Aux::String(&tag))?;

                        let mut tag = String::new();
                        let mut first = true;
                        for (gene, weight) in &info.genes_hit_reverse {
                            if !first {
                                tag.push(',')
                            }
                            first = false;
                            tag.push_str(&format!("{}={:.2}", gene, weight));
                        }
                        read.push_aux(b"XR", rust_htslib::bam::record::Aux::String(&tag))?;
                        read.push_aux(b"XP", rust_htslib::bam::record::Aux::U32(
                            info.corrected_position as u32,
                        ))?;
                    }
                }
                out_bam.write(&read)?;
            }

            ii += 1
        }
        Ok(())
    }
}

fn combine_temporary_bams(
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

pub struct TreeMatcher {
    reference_to_count_trees: HashMap<String, (OurTree, Vec<String>)>,
    reference_to_aggregation_trees: HashMap<String, (OurTree, Vec<String>)>,
}

enum ReadOut {
    Weight(Vec<(u32, f64)>, Vec<(u32, f64)>),
    Filtered(Vec<u32>, Vec<u32>),
}

impl TreeMatcher {
    fn generate_chunks(&self, bam: rust_htslib::bam::IndexedReader) -> Vec<Chunk> {
        let cg = ChunkedGenome::new(&self.reference_to_aggregation_trees, bam); // can't get the ParallelBridge to work with our lifetimes.
        cg.iter().collect()
    }

    fn hits(
        &self,
        chr: &str,
        chunk_start: u32,
        chunk_stop: u32,
        read: &rust_htslib::bam::record::Record,
        quant_is_reverse: bool,
    ) -> Result<(HashSet<String>, HashSet<String>)> {
        //I think we need the chr to get the correc ttree.
        let (tree, gene_ids) = self
            .reference_to_count_trees
            .get(chr)
            .expect("Chr not found in trees");
        let blocks = read.blocks();
        let mut gene_nos_seen_match = HashSet::<String>::new();
        let mut gene_nos_seen_reverse = HashSet::<String>::new();
        for iv in blocks.iter() {
            if (iv.1 < chunk_start)
                || iv.0 >= chunk_stop
                || ((iv.0 < chunk_start) && (iv.1 >= chunk_start))
            {
                // if this block is outside of the region
                // don't count it at all.
                // if it is on a block boundary
                // only count it for the left side.
                // which is ok, since we place the blocks to the right
                // of our intervals.
                continue;
            }
            for r in tree.find(iv.0..iv.1) {
                let entry = r.data();
                let gene_no = entry.0;
                let strand = entry.1;
                let target = match (quant_is_reverse, read.is_reverse(), strand) {
                    (false, false, Strand::Plus) => &mut gene_nos_seen_match,
                    (false, false, Strand::Minus) => &mut gene_nos_seen_reverse,
                    (false, true, Strand::Plus) => &mut gene_nos_seen_reverse,
                    (false, true, Strand::Minus) => &mut gene_nos_seen_match,

                    (true, false, Strand::Plus) => &mut gene_nos_seen_reverse,
                    (true, false, Strand::Minus) => &mut gene_nos_seen_match,
                    (true, true, Strand::Plus) => &mut gene_nos_seen_match,
                    (true, true, Strand::Minus) => &mut gene_nos_seen_reverse,
                    (_, _, Strand::Unstranded) => &mut gene_nos_seen_match,
                };
                target.insert(gene_ids[gene_no as usize].clone());
            }
        }
        Ok((gene_nos_seen_match, gene_nos_seen_reverse))
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

pub struct TagMatcher {
    tag: [u8; 2],
}

impl TagMatcher {
    fn generate_chunks(&self, bam: rust_htslib::bam::IndexedReader) -> Vec<Chunk> {
        bam.header()
            .target_names()
            .iter()
            .enumerate()
            .map(|(tid, name)| {
                let name = String::from_utf8(name.to_vec())
                    .expect("Failed to convert target name to string");
                Chunk::new(
                    name,
                    tid as u32,
                    0,
                    bam.header().target_len(tid as u32).unwrap_or(0) as u32,
                )
            })
            .collect()
    }

    fn hits(
        &self,
        _chr: &str,
        _start: u32,
        _stop: u32,
        read: &rust_htslib::bam::Record,
        _quant_is_reverse: bool,
    ) -> Result<(HashSet<String>, HashSet<String>)> {
        let mut genes_hit_correct = HashSet::new();
        let genes_hit_reverse = HashSet::new();
        if let Ok(rust_htslib::bam::record::Aux::String(value)) = read.aux(&self.tag) {
            genes_hit_correct.insert(value.to_string());
        }
        Ok((genes_hit_correct, genes_hit_reverse))
    }
}
