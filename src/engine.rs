use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
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
pub enum AnnotatedRead {
    Filtered,
    //FilteredInQuant, might want to do tthis for strategy.multi_region hits?
    NotInRegion,
    Counted(AnnotatedReadInfo),
    Duplicate,
    BarcodeNotInWhitelist,
}

#[derive(Debug)]
pub struct AnnotatedReadInfo {
    pub corrected_position: u32, // clipping corrected position
    pub genes_hit_correct: HashMap<String, f32>, //todo: optimize!
    pub genes_hit_reverse: HashMap<String, f32>, //todo: optimize!
    pub umi: Option<Vec<u8>>,    // Optional: What's it's UMI.
    pub barcode: Option<Vec<u8>>, // Optional: What's it's cell-barcode
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
    ) -> Result<(HashSet<String>, HashSet<String>)> {
        match self {
            ReadToGeneMatcher::TreeMatcher(matcher) => {
                matcher.hits(chr, chunk_start, chunk_stop, read)
            }
            ReadToGeneMatcher::TagMatcher(matcher) => {
                matcher.hits(chr, chunk_start, chunk_stop, read)
            }
        }
    }
}

pub enum Output {
    PerRegion {
        output_filename: PathBuf,
        counter: HashMap<String, (f64, f64)>,
        stat_counter: HashMap<String, usize>,
        sorted_keys: Option<Vec<String>>,
        first_column_only: bool,
        id_attribute: String,
    },
}

impl Output {
    pub fn new_per_region(
        output_filename: PathBuf,
        first_column_only: bool,
        sorted_keys: Option<Vec<String>>,
        id_attribute: String,
    ) -> Self {
        let stat_counter = [
            ("ambiguous", 0),
            ("correct", 0),
            ("reverse", 0),
            ("outside", 0),
            ("duplicate", 0),
            ("barcode_not_in_whitelist", 0),
            ("filtered", 0),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect();

        Output::PerRegion {
            output_filename,
            counter: HashMap::new(),
            stat_counter,
            sorted_keys,
            first_column_only,
            id_attribute,
        }
    }

    fn count_reads(&mut self, annotated_reads: &Vec<(AnnotatedRead, usize)>) -> Result<()> {
        match self {
            Output::PerRegion {
                counter,
                stat_counter,
                ..
            } => {
                for (read, _org_index) in annotated_reads {
                    match read {
                        AnnotatedRead::Counted(info) => {
                            match (
                                info.genes_hit_correct.is_empty(),
                                info.genes_hit_reverse.is_empty(),
                            ) {
                                (true, true) => {
                                    *stat_counter.get_mut("outside").unwrap() += 1;
                                }
                                (false, true) => *stat_counter.get_mut("correct").unwrap() += 1,
                                (true, false) => *stat_counter.get_mut("reverse").unwrap() += 1,
                                (false, false) => *stat_counter.get_mut("ambiguous").unwrap() += 1,
                            }
                            for (gene, weight) in &info.genes_hit_correct {
                                let entry = counter.entry(gene.to_string()).or_insert((0.0, 0.0));
                                entry.0 += *weight as f64; // count forward
                            }
                            for (gene, weight) in &info.genes_hit_reverse {
                                let entry = counter.entry(gene.to_string()).or_insert((0.0, 0.0));
                                entry.1 += *weight as f64; // count reverse
                            }
                        }
                        AnnotatedRead::Filtered => *stat_counter.get_mut("filtered").unwrap() += 1,
                        AnnotatedRead::NotInRegion => {
                            //these are 'accidential overfetches'
                            //*stat_counter.get_mut("outside").unwrap() += 1
                        }
                        AnnotatedRead::Duplicate => {
                            *stat_counter.get_mut("duplicate").unwrap() += 1
                        }
                        AnnotatedRead::BarcodeNotInWhitelist => {
                            *stat_counter.get_mut("barcode_not_in_whitelist").unwrap() += 1
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<()> {
        match self {
            Output::PerRegion {
                output_filename,
                counter,
                stat_counter,
                first_column_only,
                sorted_keys,
                id_attribute,
            } => {
                ex::fs::create_dir_all(&output_filename.parent().unwrap())?;
                let sorted_keys = sorted_keys.unwrap_or_else(|| {
                    let mut keys: Vec<_> = counter.keys().cloned().collect();
                    keys.sort();
                    keys
                });

                let output_file = ex::fs::File::create(&output_filename)?;
                let mut out_buffer = std::io::BufWriter::new(output_file);
                if first_column_only {
                    out_buffer
                        .write_all(format!("{}\tcount\n", id_attribute).as_bytes())
                        .context("Failed to write header to output file")?;
                    for key in sorted_keys {
                        let count = counter.get(&key).unwrap_or(&(0.0, 0.0)).0;
                        out_buffer
                            .write_all(format!("{}\t{}\n", key, count).as_bytes())
                            .context("Failed to write counts to output file")?;
                    }
                } else {
                    out_buffer
                        .write_all(
                            format!("{}\tcount_correct\tcount_reverse\n", id_attribute).as_bytes(),
                        )
                        .context("Failed to write header to output file")?;

                    for key in sorted_keys {
                        let (count_correct, count_reverse) =
                            counter.get(&key).unwrap_or(&(0.0, 0.0));
                        out_buffer
                            .write_all(
                                format!("{}\t{}\t{}\n", key, count_correct, count_reverse)
                                    .as_bytes(),
                            )
                            .context("Failed to write counts to output file")?;
                    }
                }

                Self::write_stats(&output_filename, &stat_counter)
                    .context("Failed to write stats file")?;
            }
        }
        Ok(())
    }

    fn write_stats(output_filename: &PathBuf, stat_counter: &HashMap<String, usize>) -> Result<()> {
        let stat_filename = output_filename.with_file_name(format!(
            "{}.stats.tsv",
            output_filename.file_name().unwrap().to_string_lossy()
        ));
        let output_file = ex::fs::File::create(&stat_filename)?;
        let mut out_buffer = std::io::BufWriter::new(output_file);
        out_buffer
            .write_all(b"stat\tcount\n")
            .context("Failed to write header to stats file")?;
        let sorted_keys = stat_counter.keys().sorted();
        for key in sorted_keys {
            out_buffer
                .write_all(format!("{}\t{}\n", key, stat_counter.get(key).unwrap()).as_bytes())
                .context("Failed to write stats to file")?;
        }
        Ok(())
    }
}

pub struct Engine {
    quantifier: Quantification,
    filters: Vec<crate::filters::Filter>,
    matcher: ReadToGeneMatcher,
    umi_extractor: crate::extractors::UMIExtraction,
    cell_barcode: Option<crate::barcodes::CellBarcodes>,
    output: Arc<Mutex<Output>>,
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
        cell_barcode: Option<crate::barcodes::CellBarcodes>,
        count_strategy: crate::config::Strategy,
        output: Output,
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
                count_strategy,
            }),
            filters,
            quantifier,
            umi_extractor,
            cell_barcode,
            output: Arc::new(Mutex::new(output)),
        })
    }
    pub fn from_references(
        references: Vec<(String, u64)>,
        filters: Vec<crate::filters::Filter>,
        quantifier: Quantification,
        umi_extractor: crate::extractors::UMIExtraction,
        cell_barcode: Option<crate::barcodes::CellBarcodes>,
        count_strategy: crate::config::Strategy,
        output: Output,
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
                count_strategy,
            }),
            filters,
            quantifier,
            umi_extractor,
            cell_barcode,
            output: Arc::new(Mutex::new(output)),
        })
    }

    pub fn from_bam_tag(
        tag: [u8; 2],
        filters: Vec<crate::filters::Filter>,
        quantifier: Quantification,
        umi_extractor: crate::extractors::UMIExtraction,
        cell_barcode: Option<crate::barcodes::CellBarcodes>,
        output: Output,
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
            cell_barcode,
            output: Arc::new(Mutex::new(output)),
        }
    }

    pub fn quantify_bam(
        self,
        bam_path: impl AsRef<Path>,
        index_path: Option<&Path>,
        output_path: impl AsRef<Path>,
        write_output_bam: bool,
        max_skip_len: u32,
        correct_reads_for_clipping: bool,
    ) -> Result<()> {
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
        if chunks.is_empty() {
            bail!("No chunks generated. This might be because the BAM file is empty or the matcher did not find any regions to quantify.");
        }
        //perform the counting
        let chunk_names = chunks.iter().map(|c| c.str_id()).collect::<Vec<_>>();

        let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
        let aggregated = pool.install(|| {
            let result: Vec<Result<()>> = chunks
                .into_par_iter()
                .map(|chunk| -> Result<()> {
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
                            correct_reads_for_clipping,
                        )
                        .expect("Failure in quantification");

                    //annotated_reads.retain(|read| matches!(read, AnnotatedRead::Counted(_)));
                    ////this is lifetime trouble
                    /* annotated_reads.sort_by_key(|(read, _org_index)| match read {
                        AnnotatedRead::Counted(info) => {
                            (info.corrected_position as u64, info.barcode.as_ref())
                        }
                        _ => (u32::MAX as u64 + 1, None),
                    }); */
                    annotated_reads.sort_by(|a, b| match (&a.0, &b.0) {
                        (AnnotatedRead::Counted(info_a), AnnotatedRead::Counted(info_b)) => info_a
                            .corrected_position
                            .cmp(&info_b.corrected_position)
                            .then_with(|| info_a.barcode.as_ref().cmp(&info_b.barcode.as_ref())),
                        (AnnotatedRead::Counted(_), _) => std::cmp::Ordering::Less,
                        (_, AnnotatedRead::Counted(_)) => std::cmp::Ordering::Greater,
                        _ => std::cmp::Ordering::Equal,
                    });

                    let mut last_pos = None;
                    //we can't do it without storing change_indices because we'd borrow
                    //annotated_reads twice...
                    let mut change_indices = Vec::new();
                    let mut found_filtered = false;
                    for (ii, (read, _org_index)) in annotated_reads.iter().enumerate() {
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

                    let lock = self.output.lock();
                    match lock {
                        Ok(mut output) => {
                            output
                                .count_reads(&annotated_reads)
                                .context("Failed to count reads")?;
                        }
                        Err(_) => bail!("Another thread panicked, output no longer available."),
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
                    Ok(())
                })
                .collect();
            result
        });

        if aggregated.iter().any(|r| r.is_err()) {
            let errors: Vec<_> = aggregated.into_iter().filter_map(Result::err).collect();
            bail!("Errors occurred during quantification: {:?}", errors);
        }

        let output = Arc::into_inner(self.output).context("Failed to retrieve output from arc")?;
        let output = output
            .into_inner()
            .context("Failed to unlock output mutex")?;
        output.finish()?;

        if let Some((temp_dir, output_header)) = output_bam_info {
            combine_temporary_bams(&chunk_names, temp_dir, output_prefix, output_header)?;
        }

        Ok(())
    }

    fn annotate_reads(
        &self,
        bam: &mut rust_htslib::bam::IndexedReader,
        chr: &str,
        tid: u32,
        start: u32,
        stop: u32,
        max_skip_len: u32,
        correct_reads_for_clipping: bool,
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

            let corrected_position = if correct_reads_for_clipping {
                read.corrected_pos(max_skip_len)
            } else {
                Some(read.pos())
            };

            if let Some(corrected_position) = corrected_position {
                if (
                    (corrected_position as u32) < start && start > 0
                    //reads in the first chunk, that would've started to the left are
                    //also accepted
                ) || ((corrected_position as u32) >= stop)
                {
                    //this ensures we count a read only in the chunk where it's left most
                    //pos is in.
                    res.push((AnnotatedRead::NotInRegion, ii));
                    ii += 1;
                    continue;
                }

                for f in self.filters.iter() {
                    if f.remove_read(&read) {
                        // if the read does not pass the filter, skip it
                        res.push((AnnotatedRead::Filtered, ii));
                        ii += 1;
                        continue 'outer;
                    }
                }

                let (genes_hit_correct, genes_hit_reverse) =
                    self.matcher.hits(chr, start, stop, &read)?;
                let barcode = {
                    match { self.cell_barcode.as_ref() } {
                        Some(cb) => {
                            match cb.correct(&cb.extract(&read).expect("barcode extraction failed"))
                            {
                                Some(corrected_barcode) => Some(corrected_barcode),
                                None => {
                                    res.push((AnnotatedRead::BarcodeNotInWhitelist, ii));
                                    ii += 1;
                                    continue;
                                }
                            }
                        }
                        None => None,
                    }
                };

                let info = AnnotatedReadInfo {
                    corrected_position: corrected_position as u32,
                    genes_hit_correct: genes_hit_correct.into_iter().map(|id| (id, 1.0)).collect(),
                    genes_hit_reverse: genes_hit_reverse.into_iter().map(|id| (id, 1.0)).collect(),
                    umi: self.umi_extractor.extract(&read),
                    barcode: barcode,
                    mapping_priority: (
                        read.no_of_alignments().try_into().unwrap_or(255),
                        read.mapq(),
                    ),
                };
                res.push((AnnotatedRead::Counted(info), ii));
                ii += 1;
            } else {
                panic!("Did not expect to see an unaligned read in a chunk: {:?}, pos: {}, cigar: {} - corrected to {:?}", 
                    std::str::from_utf8(read.qname()).unwrap(), read.pos(), read.cigar(), read.corrected_pos(max_skip_len));
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
                    AnnotatedRead::Filtered => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(1))?;
                    }
                    /* AnnotatedRead::FilteredInQuant => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(2))?;
                    } */
                    AnnotatedRead::Duplicate => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(3))?;
                    }
                    AnnotatedRead::BarcodeNotInWhitelist => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(4))?;
                    }
                    AnnotatedRead::Counted(info) => {
                        //we have a read that was annotated
                        //write it to the output bam
                        let mut tag = String::new();
                        let mut first = true;

                        for gene in info.genes_hit_correct.keys().sorted() {
                            let weight = *info.genes_hit_correct.get(gene).unwrap();
                            if !first {
                                tag.push(',')
                            }
                            first = false;
                            tag.push_str(&format!("{}={:.2}", gene, weight));
                        }
                        read.replace_aux(b"XQ", rust_htslib::bam::record::Aux::String(&tag))?;

                        let mut tag = String::new();
                        let mut first = true;
                        for gene in info.genes_hit_reverse.keys().sorted() {
                            let weight = *info.genes_hit_reverse.get(gene).unwrap();
                            if !first {
                                tag.push(',')
                            }
                            first = false;
                            tag.push_str(&format!("{}={:.2}", gene, weight));
                        }
                        read.replace_aux(b"XR", rust_htslib::bam::record::Aux::String(&tag))?;
                        read.replace_aux(
                            b"XP",
                            rust_htslib::bam::record::Aux::U32(info.corrected_position as u32),
                        )?;

                        if let Some(cell_barcode) = info.barcode.as_ref() {
                            read.replace_aux(
                                b"CB",
                                rust_htslib::bam::record::Aux::String(
                                    std::str::from_utf8(cell_barcode).expect("barcode wasn't utf8"),
                                ),
                            )?;
                        }
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
    count_strategy: crate::config::Strategy,
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
    ) -> Result<(HashSet<String>, HashSet<String>)> {
        use crate::config::MatchDirection;
        let (tree, gene_ids) = self
            .reference_to_count_trees
            .get(chr)
            .expect("Chr not found in trees");
        let blocks = read.blocks();
        let mut gene_nos_seen_match = HashMap::<String, HashSet<u32>>::new();
        let mut gene_nos_seen_reverse = HashMap::<String, HashSet<u32>>::new();
        let mut bases_aligned = 0u32;
        for iv in blocks.iter() {
            bases_aligned += iv.1 - iv.0;
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
                let found_interval = r.interval();
                let overlap_range = match self.count_strategy.overlap {
                    crate::config::OverlapMode::Union => HashSet::new(),
                    _ => (iv.0.max(found_interval.start)..iv.1.min(found_interval.end)).collect(),
                };

                let gene_no = entry.0;
                let region_strand = entry.1;
                let target = match (
                    &self.count_strategy.direction,
                    read.is_reverse(),
                    region_strand,
                ) {
                    (MatchDirection::Forward, false, Strand::Plus) => &mut gene_nos_seen_match,
                    (MatchDirection::Forward, false, Strand::Minus) => &mut gene_nos_seen_reverse,
                    (MatchDirection::Forward, true, Strand::Plus) => &mut gene_nos_seen_reverse,
                    (MatchDirection::Forward, true, Strand::Minus) => &mut gene_nos_seen_match,
                    (MatchDirection::Forward, _, Strand::Unstranded) => &mut gene_nos_seen_match,

                    (MatchDirection::Reverse, false, Strand::Plus) => &mut gene_nos_seen_reverse,
                    (MatchDirection::Reverse, false, Strand::Minus) => &mut gene_nos_seen_match,
                    (MatchDirection::Reverse, true, Strand::Plus) => &mut gene_nos_seen_match,
                    (MatchDirection::Reverse, true, Strand::Minus) => &mut gene_nos_seen_reverse,
                    (MatchDirection::Reverse, _, Strand::Unstranded) => &mut gene_nos_seen_match,
                    (MatchDirection::Ignore, _, _) => &mut gene_nos_seen_match,
                };
                target
                    .entry(gene_ids[gene_no as usize].clone())
                    .and_modify(|e| e.extend(&overlap_range))
                    .or_insert(overlap_range);
            }
        }

        for gg in [&mut gene_nos_seen_match, &mut gene_nos_seen_reverse] {
            match self.count_strategy.overlap {
                crate::config::OverlapMode::Union => {
                    // do nothing, we keep them all
                }
                crate::config::OverlapMode::IntersectionStrict => {
                    //only keep those that are fully contained in the region
                    if std::str::from_utf8(read.qname()).unwrap()
                        == "HWI-C00113:73:HVM2VBCXX:1:1209:15870:92336"
                    {
                        dbg!(std::str::from_utf8(read.qname()).unwrap(), &gg);
                    }
                    gg.retain(|_, v| v.len() == bases_aligned as usize);
                }
                crate::config::OverlapMode::IntersectionNonEmpty => {
                    let any_fully_contained =
                        gg.values().any(|v| v.len() == bases_aligned as usize);
                    if any_fully_contained {
                        //only keep those that are fully contained in the region
                        gg.retain(|_, v| v.len() == bases_aligned as usize);
                    } else {
                        // keep all.
                    }
                }
            }
            match self.count_strategy.multi_region {
                crate::config::MultiRegionHandling::Drop => {
                    if gg.len() > 1 {
                        // if there are multiple genes, drop them
                        gg.clear();
                    }
                }
                crate::config::MultiRegionHandling::CountBoth => {
                    //do nothing.
                }
            }
        }

        Ok((
            gene_nos_seen_match.into_keys().collect(),
            gene_nos_seen_reverse.into_keys().collect(),
        ))
    }
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
    ) -> Result<(HashSet<String>, HashSet<String>)> {
        let mut genes_hit_correct = HashSet::new();
        let genes_hit_reverse = HashSet::new();
        if let Ok(rust_htslib::bam::record::Aux::String(value)) = read.aux(&self.tag) {
            genes_hit_correct.insert(value.to_string());
        }
        Ok((genes_hit_correct, genes_hit_reverse))
    }
}
