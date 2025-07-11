use ex::Wrapper;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use string_interner::symbol::SymbolU32;
use string_interner::StringInterner;
mod chunked_genome;

use crate::bam_ext::BamRecordExtensions;
use crate::config::MatchDirection;
use crate::deduplication::{AcceptReadResult, DeduplicationStrategy};
use crate::extractors::UMIExtractor;
use crate::filters::ReadFilter;
use crate::gtf::Strand;
use anyhow::{bail, Context, Result};
use bio::data_structures::interval_tree::IntervalTree;
use chunked_genome::{Chunk, ChunkedGenome};
use itertools::{izip, Itertools};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_htslib::bam::{self, Read as ReadTrait};

use crate::gtf::GTFEntrys;
pub type OurTree = IntervalTree<u32, (u32, Strand)>;

pub type OurInterner = StringInterner<string_interner::backend::StringBackend>;

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
        let by_chr_interval = intervals.entry(*seq_name_cat_id).or_default();
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
    Counted(Box<AnnotatedReadInfo>), // otherwise this is as large as AnnotatedReadInfo for each
    // entry.
    Duplicate,
    NoBarcode,
    NoUMI,
    BarcodeNotInWhitelist(Vec<u8>),
}

#[derive(Debug)]
pub struct Hits {
    correct: Vec<string_interner::symbol::SymbolU32>,
    reverse: Vec<string_interner::symbol::SymbolU32>,
}

#[derive(Debug)]
pub struct AnnotatedReadInfo {
    pub corrected_position: i32, // clipping corrected position. Samspec is limited to 0..2^31-1,
    // so i32 is safe, and it allows us to have negative corrected positions
    pub hits: Hits,
    pub umi: Option<Vec<u8>>,     // Optional: What's it's UMI. 24 bytes
    pub barcode: Option<Vec<u8>>, // Optional: What's it's cell-barcode 24 bytes
    pub mapping_priority: (u8, u8),
    pub reverse: bool,
}

pub enum ReadToGeneMatcher {
    TreeMatcher(TreeMatcher),
    TagMatcher(TagMatcher),
    ReferenceMatcher(ReferenceMatcher),
}

impl ReadToGeneMatcher {
    fn generate_chunks(&self, bam: rust_htslib::bam::IndexedReader) -> Result<Vec<Chunk>> {
        match self {
            ReadToGeneMatcher::TreeMatcher(matcher) => matcher.generate_chunks(bam),
            ReadToGeneMatcher::TagMatcher(matcher) => matcher.generate_chunks(bam),
            ReadToGeneMatcher::ReferenceMatcher(matcher) => matcher.generate_chunks(bam),
        }
    }

    fn hits(
        &self,
        chunk: &Chunk,
        read: &rust_htslib::bam::record::Record,
        interner: &mut OurInterner,
    ) -> Result<(
        Vec<string_interner::symbol::SymbolU32>,
        Vec<string_interner::symbol::SymbolU32>,
    )> {
        match self {
            ReadToGeneMatcher::TreeMatcher(matcher) => matcher.hits(chunk, read, interner),
            ReadToGeneMatcher::TagMatcher(matcher) => matcher.hits(chunk, read, interner),
            ReadToGeneMatcher::ReferenceMatcher(matcher) => matcher.hits(chunk, read, interner),
        }
    }
}

pub enum CounterPerChunk {
    PerRegion {
        counter: HashMap<string_interner::symbol::SymbolU32, (usize, usize)>,
        stat_counter: HashMap<String, usize>,
    },
    SingleCell {
        stat_counter: HashMap<String, usize>,
        counter: HashMap<(string_interner::symbol::SymbolU32, Vec<u8>), usize>,
    },
}

impl CounterPerChunk {
    fn count_reads(&mut self, annotated_reads: &Vec<(AnnotatedRead, usize)>) -> Result<()> {
        match self {
            CounterPerChunk::PerRegion {
                counter,
                stat_counter,
                ..
            } => {
                for (read, _org_index) in annotated_reads {
                    let count_as = match read {
                        AnnotatedRead::Counted(info) => {
                            let hits = &info.hits;
                            for gene in &hits.correct {
                                let entry = counter.entry(*gene).or_insert((0, 0));
                                entry.0 = entry.0.saturating_add(1)
                            }
                            for gene in &hits.reverse {
                                let entry = counter.entry(*gene).or_insert((0, 0));
                                entry.1 = entry.1.saturating_add(1)
                            }
                            match (hits.correct.is_empty(), hits.reverse.is_empty()) {
                                (true, true) => "outside",
                                (false, true) => "correct",
                                (true, false) => "reverse",
                                (false, false) => "ambiguous",
                            }
                        }
                        AnnotatedRead::Filtered => "filtered",
                        AnnotatedRead::NotInRegion => {
                            //these are 'accidential overfetches'
                            "overfetch"
                        }
                        AnnotatedRead::Duplicate => "duplicate",
                        AnnotatedRead::NoBarcode => "no_barcode",

                        AnnotatedRead::NoUMI => "no_umi",
                        AnnotatedRead::BarcodeNotInWhitelist(_) => "barcode_not_in_whitelist",
                    };
                    if count_as != "overfetch" {
                        //todo: preinsert values
                        *stat_counter.entry(count_as.to_string()).or_default() += 1;
                    }
                }
            }
            CounterPerChunk::SingleCell {
                stat_counter,
                counter,
            } => {
                //we want a consistent output order
                for (read, _org_index) in annotated_reads {
                    let count_as = match read {
                        AnnotatedRead::Counted(info) => {
                            let barcode = info.barcode.as_ref().expect("no barcode? bug");
                            let hits = &info.hits;
                            for gene in &hits.correct {
                                counter
                                    .entry((*gene, barcode.clone()))
                                    .and_modify(|e| *e += 1)
                                    .or_insert(1);
                            }
                            match (hits.correct.is_empty(), hits.reverse.is_empty()) {
                                (true, true) => "outside",
                                (false, true) => "correct",
                                (true, false) => "reverse",
                                (false, false) => "ambiguous",
                            }
                        }
                        AnnotatedRead::Filtered => "filtered",
                        AnnotatedRead::NotInRegion =>
                        //these are 'accidential overfetches'
                        {
                            "overfetch"
                        }

                        AnnotatedRead::Duplicate => "duplicate",

                        AnnotatedRead::NoBarcode => "no_barcode",

                        AnnotatedRead::NoUMI => "no_umi",
                        AnnotatedRead::BarcodeNotInWhitelist(_) => "barcode_not_in_whitelist",
                    };
                    if count_as != "overfetch" {
                        //todo: preinsert values
                        *stat_counter.entry(count_as.to_string()).or_default() += 1;
                    }
                }
            }
        }
        Ok(())
    }
}

pub enum Output {
    PerRegion {
        output_filename: PathBuf,
        counter: HashMap<String, (usize, usize)>,
        stat_counter: HashMap<String, usize>,
        sorted_keys: Option<Vec<String>>,
        first_column_only: bool,
        id_attribute: String,
    },
    SingleCell {
        output_prefix: PathBuf,
        stat_counter: HashMap<String, usize>,
        features: HashMap<String, usize>,
        barcodes: HashSet<Vec<u8>>,
        matrix_temp_dir: PathBuf,
        entry_count: usize,
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
            ("no_barcode", 0),
            ("no_umi", 0),
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

    pub fn new_singlecell(
        output_prefix: PathBuf,
        sorted_keys: Option<Vec<String>>,
    ) -> Result<Self> {
        let matrix_temp_dir = output_prefix.join("matrix.mtx.temp");
        if matrix_temp_dir.exists() {
            //remove the whole tree
            ex::fs::remove_dir_all(&matrix_temp_dir)
                .context("Failed to remove existing temp directory")?;
        }
        ex::fs::create_dir_all(&matrix_temp_dir)
            .context("failed to create output directory (or the nested temp dir)")?;

        let mut features = HashMap::new();
        if let Some(sorted_keys) = sorted_keys.as_ref() {
            for key in sorted_keys {
                features.insert(key.clone(), features.len());
            }
        }
        let stat_counter = [
            ("ambiguous", 0),
            ("correct", 0),
            ("reverse", 0),
            ("outside", 0),
            ("duplicate", 0),
            ("no_barcode", 0),
            ("no_umi", 0),
            ("barcode_not_in_whitelist", 0),
            ("filtered", 0),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect();

        Ok(Output::SingleCell {
            output_prefix,
            stat_counter,
            features,
            barcodes: HashSet::new(),
            matrix_temp_dir,
            entry_count: 0,
        })
    }

    pub fn per_chunk(&self) -> CounterPerChunk {
        match self {
            Output::PerRegion { .. } => CounterPerChunk::PerRegion {
                counter: HashMap::new(),
                stat_counter: HashMap::new(),
            },
            Output::SingleCell { .. } => CounterPerChunk::SingleCell {
                counter: HashMap::new(),
                stat_counter: HashMap::new(),
            },
        }
    }

    //the final, we're done with this chunk, push some output
    fn count_reads(
        &mut self,
        chunk: &Chunk,
        interner: &OurInterner,
        output_catcher: CounterPerChunk,
    ) -> Result<()> {
        match self {
            Output::PerRegion {
                counter,
                stat_counter,
                ..
            } => {
                if let CounterPerChunk::PerRegion {
                    counter: incoming_counter,
                    stat_counter: incoming_stat_counter,
                } = output_catcher
                {
                    for (k, v) in incoming_counter.iter() {
                        let entry = counter
                            .entry(interner.resolve(*k).unwrap().to_string())
                            .or_insert((0, 0));
                        entry.0 += v.0;
                        entry.1 += v.1;
                    }
                    for (k, v) in incoming_stat_counter {
                        *stat_counter.entry(k.clone()).or_default() += v;
                    }
                } else {
                    unreachable!();
                }
            }
            Output::SingleCell {
                matrix_temp_dir,
                features,
                barcodes,
                stat_counter,
                entry_count,
                ..
            } => {
                //we want a consistent output order

                if let CounterPerChunk::SingleCell {
                    counter: incoming_counter,
                    stat_counter: incoming_stat_counter,
                } = output_catcher
                {
                    for (k, v) in incoming_stat_counter {
                        *stat_counter.entry(k.clone()).or_default() += v;
                    }

                    let mut matrix_handle = BufWriter::new(
                        ex::fs::File::create(matrix_temp_dir.join(chunk.str_id()))?.into_inner(),
                    );

                    //now these genes are fully measured, write them out
                    for ((feature_ref, barcode), value) in incoming_counter.into_iter() {
                        let feature_str = interner.resolve(feature_ref).unwrap();
                        let features_len = features.len();
                        let feature_idx = match features.entry(feature_str.to_string()) {
                            std::collections::hash_map::Entry::Occupied(e) => *e.get(),
                            std::collections::hash_map::Entry::Vacant(e) => {
                                let new_index = features_len;
                                e.insert(new_index);
                                new_index
                            }
                        };
                        barcodes.insert(barcode.clone());
                        *entry_count += 1;
                        matrix_handle
                            .write_all(
                                format!(
                                    "{} {} {}\n",
                                    feature_idx + 1, // MatrixMarket is 1-based
                                    std::str::from_utf8(&barcode).unwrap(),
                                    value
                                )
                                .as_bytes(),
                            )
                            .context("Failed to write to matrix file")?;
                    }
                }
            }
        }
        Ok(())
    }

    fn finish(self, chunk_names: &[String]) -> Result<()> {
        measure_time::info_time!("Preparing final output");
        match self {
            Output::PerRegion {
                output_filename,
                counter,
                stat_counter,
                first_column_only,
                sorted_keys,
                id_attribute,
            } => {
                ex::fs::create_dir_all(output_filename.parent().unwrap())?;
                let sorted_keys = sorted_keys.unwrap_or_else(|| {
                    let mut keys: Vec<_> = counter.keys().map(|x| x.to_string()).collect();
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
                        let count = counter.get(&key).unwrap_or(&(0, 0)).0;
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
                        let (count_correct, count_reverse) = counter.get(&key).unwrap_or(&(0, 0));
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
            Output::SingleCell {
                output_prefix,
                features,
                barcodes,
                matrix_temp_dir,
                entry_count,
                stat_counter,
            } => {
                let features_filename = output_prefix.join("features.tsv");
                let feature_len = features.len();
                let sorted_features = {
                    let mut temp = features
                        .into_iter()
                        .map(|(feature, index)| (index, feature))
                        .collect::<Vec<_>>();
                    temp.sort();
                    temp
                };
                let mut feature_file =
                    BufWriter::new(ex::fs::File::create(&features_filename)?.into_inner());
                for (_, feature) in sorted_features {
                    feature_file
                        .write_all(format!("{}\n", feature).as_bytes())
                        .context("Failed to write features to file")?;
                }

                let barcodes_filename = output_prefix.join("barcodes.tsv");
                let barcode_len = barcodes.len();
                let (barcodes, barcode_to_index) = {
                    let mut temp: Vec<_> = barcodes.into_iter().collect();
                    temp.sort();
                    let lookup: HashMap<Vec<u8>, usize> = temp
                        .iter()
                        .enumerate()
                        .map(|(i, b)| (b.clone(), i + 1))
                        .collect();

                    (temp, lookup)
                };
                let mut barcode_file =
                    BufWriter::new(ex::fs::File::create(&barcodes_filename)?.into_inner());
                for barcode in barcodes {
                    barcode_file
                        .write_all(format!("{}\n", String::from_utf8_lossy(&barcode)).as_bytes())
                        .context("Failed to write barcodes to file")?;
                }
                //matrix file has a header with nrow, ncols, nentries, so we need to push it into a
                //new file.
                let matrix_filename = output_prefix.join("matrix.mtx");

                let mut matrix_out = ex::fs::File::create(&matrix_filename)?;
                matrix_out.write_all(
                    "%%MatrixMarket matrix coordinate integer general\n%\n".as_bytes(),
                )?;
                matrix_out
                    .write_all(
                        format!("{} {} {}\n", feature_len, barcode_len, entry_count,).as_bytes(),
                    )
                    .context("Failed to write header to matrix file")?;

                let sorted_chunk_names: Vec<_> = chunk_names.iter().sorted().collect();

                for chunk_str_id in sorted_chunk_names {
                    let temp_filename = matrix_temp_dir.join(chunk_str_id);
                    let temp_handle = BufReader::new(
                        ex::fs::File::open(&temp_filename)
                            .context("Failed to open temporary matrix file")?
                            .into_inner(),
                    );
                    for line in temp_handle.lines() {
                        let line =
                            line.context("Failed to read line from temporary matrix file")?;
                        let mut parts = line.split_whitespace();
                        let feature_idx: usize = parts
                            .next()
                            .context("Missing feature index in line")?
                            .parse()
                            .context("Failed to parse feature index")?;
                        let barcode = parts
                            .next()
                            .context("Missing barcode in line")?
                            .as_bytes()
                            .to_vec();
                        let value: f64 = parts
                            .next()
                            .context("Missing value in line")?
                            .parse()
                            .context("Failed to parse value")?;

                        if let Some(&barcode_index) = barcode_to_index.get(&barcode) {
                            matrix_out.write_all(
                                format!(
                                    "{} {} {}\n",
                                    feature_idx, // we already added +1
                                    barcode_index,
                                    value
                                )
                                .as_bytes(),
                            )?;
                        }
                    }
                }

                if matrix_temp_dir.exists() {
                    ex::fs::remove_dir_all(&matrix_temp_dir)
                        .context("Failed to remove temporary matrix directory")?;
                } else {
                    bail!(
                        "Temporary matrix dir did not exist, what did we just copy?: {}",
                        matrix_temp_dir.display()
                    );
                }

                Self::write_stats(&matrix_filename, &stat_counter)
                    .context("Failed to write stats file")?;
            }
        }
        Ok(())
    }

    fn write_stats(output_filename: &Path, stat_counter: &HashMap<String, usize>) -> Result<()> {
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

struct PerPosition {
    reads_forward: Vec<(AnnotatedRead, usize)>,
    reads_reverse: Vec<(AnnotatedRead, usize)>,
    dedup_storage: crate::deduplication::DedupPerBucket,
}

struct OutputBamInfo {
    output_bam_path: PathBuf,
    header: rust_htslib::bam::Header,
}

#[derive(Debug, Copy, Clone)]
enum Bucket {
    PerPosition,
    PerReference(i32),
}

pub struct Engine {
    dedup_strategy: DeduplicationStrategy,
    filters: Vec<crate::filters::Filter>,
    matcher: ReadToGeneMatcher,
    umi_extractor: Option<crate::extractors::UMIExtraction>,
    cell_barcode: Option<crate::barcodes::CellBarcodes>,
    output: Arc<Mutex<Output>>,
}

impl Engine {
    #[allow(clippy::too_many_arguments)]
    pub fn from_gtf(
        mut gtf_entries: HashMap<String, GTFEntrys>,
        entry_kind: &str,
        entry_id_attribute: &str,
        aggregation_id_attribute: &str,
        filters: Vec<crate::filters::Filter>,
        dedup_strategy: DeduplicationStrategy,
        umi_extractor: Option<crate::extractors::UMIExtraction>,
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
            dedup_strategy,
            umi_extractor,
            cell_barcode,
            output: Arc::new(Mutex::new(output)),
        })
    }
    pub fn from_references(
        references: Vec<(String, i32, u64)>,
        filters: Vec<crate::filters::Filter>,
        dedup_strategy: DeduplicationStrategy,
        umi_extractor: Option<crate::extractors::UMIExtraction>,
        cell_barcode: Option<crate::barcodes::CellBarcodes>,
        count_strategy: crate::config::Strategy,
        output: Output,
    ) -> Result<Self> {
        let tids_to_include: Vec<i32> = references
            .iter()
            .map(|(_name, tid, _length)| *tid)
            .collect();
        Ok(Engine {
            matcher: ReadToGeneMatcher::ReferenceMatcher(ReferenceMatcher {
                tids_to_include,
                direction: count_strategy.direction,
            }),
            filters,
            dedup_strategy,
            umi_extractor,
            cell_barcode,
            output: Arc::new(Mutex::new(output)),
        })
    }

    pub fn from_bam_tag(
        tag: [u8; 2],
        filters: Vec<crate::filters::Filter>,
        dedup_strategy: DeduplicationStrategy,
        umi_extractor: Option<crate::extractors::UMIExtraction>,
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
            dedup_strategy,
            umi_extractor,
            cell_barcode,
            output: Arc::new(Mutex::new(output)),
        }
    }

    pub fn quantify_bam(
        mut self,
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

        let bam_header = bam.header();
        for filter in self.filters.iter_mut() {
            filter.init(bam_header).context("filter init failed")?;
        }

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

                Some(OutputBamInfo {
                    output_bam_path: od,
                    header,
                })
            }

            false => None,
        };
        let chunks = {
            let mut chunks = self.matcher.generate_chunks(bam)?;
            if chunks.is_empty() {
                bail!("No chunks generated. This might be because the BAM file is empty or the matcher did not find any regions to quantify.");
            }
            for f in self.filters.iter() {
                if let crate::filters::Filter::Reference(reference_filter) = f {
                    match reference_filter.action {
                        crate::filters::KeepOrRemove::Keep => {
                            chunks.retain(|c| reference_filter.references.contains(&c.chr));
                        }
                        crate::filters::KeepOrRemove::Remove => {
                            chunks.retain(|c| !reference_filter.references.contains(&c.chr));
                        }
                    }
                }
            }
            chunks
        };

        //perform the counting
        let chunk_names = chunks.iter().map(|c| c.str_id()).collect::<Vec<_>>();

        let pool = rayon::ThreadPoolBuilder::new().build().unwrap();

        let max_skip_len = if correct_reads_for_clipping {
            assert!(max_skip_len < i32::MAX as u32);
            max_skip_len
        } else {
            0u32
        };
        let aggregated = pool.install(|| {
            let result: Vec<Result<()>> = chunks
                .into_par_iter()
                .map(|chunk| -> Result<()> {
                    {
                        let mut bam =
                            crate::io::open_indexed_bam(bam_filename, index_filename).unwrap();
                        // Within one chunk, the reads we see will fit the same genes,
                        // over and over,
                        // so interning the strings is a good memory saving measure.
                        let mut interner = StringInterner::default();
                        let mut read_catcher: BTreeMap<i32, PerPosition> = BTreeMap::new();
                        let mut output_catcher = self
                            .output
                            .lock()
                            .expect("Another thread panicked, output no longer available.")
                            .per_chunk();

                        let mut idx_to_annotation_decision =
                            output_bam_info.as_ref().map(|_| HashMap::new());

                        let mut current_pos = 0i32;

                        bam.fetch((
                            chunk.tid,
                            chunk.start as u64,
                            (chunk.stop + max_skip_len) as u64, // this is a correctness issue
                                                                // We'll read unnecessary reads here,
                                                                // just to discard them because their corrected position is beyond the region
                                                                // but there might be a read there that is corrected to
                                                                // be within.
                        ))?;
                        let mut read = bam::Record::new();
                        let mut orig_index = 0;
                        let mut debug_processed_ids: HashSet<i32> = HashSet::new();

                        let bucket_mode = match self.dedup_strategy.bucket {
                            crate::deduplication::DeduplicationBucket::PerPosition => {
                                Bucket::PerPosition
                            }
                            crate::deduplication::DeduplicationBucket::PerReference => {
                                Bucket::PerReference((chunk.stop + max_skip_len) as i32)
                            }
                        };

                        while let Some(next_pos) = self.read_until_next_pos(
                            &mut bam,
                            &mut read,
                            &mut orig_index,
                            &chunk,
                            max_skip_len,
                            &mut interner,
                            current_pos,
                            &mut read_catcher,
                            bucket_mode,
                        )? {
                            let remaining_positions =
                                read_catcher.split_off(&(current_pos - max_skip_len as i32));
                            let done_positions = read_catcher;
                            read_catcher = remaining_positions;
                            for (done_pos, block) in done_positions.into_iter() {
                                if debug_processed_ids.contains(&done_pos) {
                                    panic!("We are processing the same position twice.Bug");
                                }
                                debug_processed_ids.insert(done_pos);
                                self.capture_read_block(
                                    block,
                                    &mut output_catcher,
                                    idx_to_annotation_decision.as_mut(),
                                )?;
                            }
                            current_pos = next_pos;
                        }
                        //capture eventual remaining blocks
                        for (done_pos, block) in read_catcher.into_iter() {
                            if debug_processed_ids.contains(&done_pos) {
                                panic!("We are processing the same position twice.Bug");
                            }
                            for read in &block.reads_forward {
                                if let AnnotatedRead::Counted(info) = &read.0 {
                                    assert_eq!(info.corrected_position, done_pos);
                                }
                            }

                            for read in &block.reads_reverse {
                                if let AnnotatedRead::Counted(info) = &read.0 {
                                    assert_eq!(info.corrected_position, done_pos);
                                }
                            }

                            self.capture_read_block(
                                block,
                                &mut output_catcher,
                                idx_to_annotation_decision.as_mut(),
                            )?;
                        }

                        self.output
                            .lock()
                            .expect("Another thread panicked, output no longer available.")
                            .count_reads(&chunk, &interner, output_catcher)
                            .context("Failed to count reads")?;

                        if let Some(OutputBamInfo {
                            output_bam_path,
                            header,
                        }) = output_bam_info.as_ref()
                        {
                            Self::write_annotated_reads(
                                &mut bam,
                                &chunk,
                                idx_to_annotation_decision.take().unwrap(),
                                output_bam_path,
                                header,
                                max_skip_len,
                                &interner,
                            )
                            .context("Failed to write output bam")?;
                        }
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
        output.finish(&chunk_names)?;

        if let Some(OutputBamInfo {
            output_bam_path,
            header,
        }) = output_bam_info
        {
            combine_temporary_bams(&chunk_names, output_bam_path, output_prefix, header)?;
            /* println!(
                "Output written to: {}",
                output_prefix.join("annotated.bam").display()
            ); */
        }

        Ok(())
    }

    fn capture_read_block(
        &self,
        read_block: PerPosition,
        output_cacher: &mut CounterPerChunk,
        mut idx_to_annotated: Option<&mut HashMap<usize, AnnotatedRead>>,
    ) -> Result<()> {
        for block in [read_block.reads_forward, read_block.reads_reverse] {
            output_cacher
                .count_reads(&block)
                .context("Failed to count reads in read block")?;
            if let Some(idx_to_annotated) = idx_to_annotated.as_mut() {
                for (read, org_index) in block {
                    idx_to_annotated.insert(org_index, read);
                }
            }
        }

        /* // we still need to establish the barcode sorting, even if
        // correct_reads_for_clipping is false
                } */

        /*{ measure_time::info_time!(
                "{}:{} - outputing reads",
                &chunk.chr,
                chunk.start
            );

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
                self.dedup_strategy
                    .dedup(&mut annotated_reads[*start..*stop])
                    .context("weighting failed")?;
            }

            let lock = self.output.lock();
            match lock {
                Ok(mut output) => {
                    output
                        .count_reads(&annotated_reads, &chunk, &interner)
                        .context("Failed to count reads")?;
                }
                Err(_) => {
                    bail!("Another thread panicked, output no longer available.")
                }
            }
        } */
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn read_until_next_pos(
        &self,
        bam: &mut rust_htslib::bam::IndexedReader,
        read: &mut bam::Record,
        org_index: &mut usize,
        chunk: &Chunk,
        max_skip_len: u32,
        interner: &mut OurInterner,
        current_pos: i32,
        read_catcher: &mut BTreeMap<i32, PerPosition>,
        bucket_mode: Bucket,
    ) -> Result<Option<i32>> {
        let mut last_read_pos: Option<i32> = None;
        'outer: loop {
            if let Some(last_read_pos) = last_read_pos {
                if last_read_pos > current_pos {
                    //we are done with this chunk
                    return Ok(Some(last_read_pos));
                }
            }

            match bam.read(read) {
                Some(Ok(result)) => result,
                Some(Err(e)) => {
                    bail!(e);
                }
                None => return Ok(None),
            };

            let (corrected_position, position_for_bounds_check) = match bucket_mode {
                Bucket::PerPosition => {
                    if max_skip_len > 0 {
                        let rp = read
                            .corrected_pos(max_skip_len)
                            .expect("unaligned read found?");
                        (rp, rp)
                    } else {
                        (read.pos() as i32, read.pos() as i32) //this is safe. it's an aligned read, must be <2^31
                    }
                }
                Bucket::PerReference(l) => (l, read.pos() as i32),
            };

            last_read_pos = Some(
                read.pos()
                    .try_into()
                    .expect("sam read position for aligned reads should always be <=2^31-1"),
            ); //we read based on stored read position. Once that's
               //reached x+max_skip len the outer frame can process those up to x...

            let output = read_catcher
                .entry(corrected_position)
                .or_insert_with(|| PerPosition {
                    reads_forward: Vec::new(),
                    reads_reverse: Vec::new(),
                    dedup_storage: self.dedup_strategy.new_bucket(),
                });
            let res = if read.is_reverse() {
                &mut output.reads_reverse
            } else {
                &mut output.reads_forward
            };
            if (chunk.start > 0 && position_for_bounds_check < chunk.start as i32)
                    //reads in the first chunk, that would've started to the left are
                    //also accepted
                 || (position_for_bounds_check >= chunk.stop as i32)
            {
                //this ensures we count a read only in the chunk where it's left most
                //pos is in.
                dbg!(
                    "not in region",
                    std::str::from_utf8(read.qname()).unwrap(),
                    corrected_position,
                    chunk
                );
                res.push((AnnotatedRead::NotInRegion, *org_index));
                *org_index += 1;
                continue;
            }

            for f in self.filters.iter() {
                if f.remove_read(read) {
                    // if the read does not pass the filter, skip it
                    res.push((AnnotatedRead::Filtered, *org_index));
                    *org_index += 1;
                    continue 'outer;
                }
            }

            let barcode = {
                match self.cell_barcode.as_ref() {
                    Some(cb) => {
                        let bc = cb.extract(read).context("barcode extraction failed")?; // an error
                        match bc {
                            Some(uncorrected) => {
                                // if we have a barcode, correct it
                                let corrected_barcode = cb.correct(&uncorrected);
                                match corrected_barcode {
                                    Some(bc) => Some(bc),
                                    None => {
                                        res.push((
                                            AnnotatedRead::BarcodeNotInWhitelist(
                                                uncorrected.clone(),
                                            ),
                                            *org_index,
                                        ));
                                        *org_index += 1;
                                        continue;
                                    }
                                }
                            }
                            None => {
                                res.push((AnnotatedRead::NoBarcode, *org_index));
                                *org_index += 1;
                                continue;
                            }
                        }
                    }
                    None => None,
                }
            };
            let umi: Option<Vec<u8>> = {
                match self.umi_extractor.as_ref() {
                    Some(x) => match x.extract(read)? {
                        Some(umi) => Some(umi),
                        None => {
                            res.push((AnnotatedRead::NoUMI, *org_index));
                            *org_index += 1;
                            continue;
                        }
                    },
                    None => None,
                }
            };
            let (genes_hit_correct, genes_hit_reverse) =
                self.matcher.hits(chunk, read, interner)?;

            let do_accept: AcceptReadResult =
                output
                    .dedup_storage
                    .accept_read(read, res.len(), umi.as_ref(), barcode.as_ref());
            match do_accept {
                AcceptReadResult::Duplicated => {
                    res.push((AnnotatedRead::Duplicate, *org_index));
                    *org_index += 1;
                }
                AcceptReadResult::New => {
                    let info = AnnotatedReadInfo {
                        corrected_position: corrected_position,
                        reverse: read.is_reverse(),
                        hits: Hits {
                            correct: genes_hit_correct.into_iter().collect(),
                            reverse: genes_hit_reverse.into_iter().collect(),
                        },
                        umi,
                        barcode,
                        mapping_priority: (
                            read.no_of_alignments().try_into().unwrap_or(255),
                            read.mapq(),
                        ),
                    };
                    res.push((AnnotatedRead::Counted(Box::new(info)), *org_index));
                    *org_index += 1;
                }
                AcceptReadResult::DuplicateButPrefered(idx_to_set_duplicate) => {
                    let (_old, old_org_index) = &res[idx_to_set_duplicate];
                    res[idx_to_set_duplicate] = (AnnotatedRead::Duplicate, *old_org_index);

                    let info = AnnotatedReadInfo {
                        corrected_position: corrected_position,
                        reverse: read.is_reverse(),
                        hits: Hits {
                            correct: genes_hit_correct,
                            reverse: genes_hit_reverse,
                        },
                        umi,
                        barcode,
                        mapping_priority: (
                            read.no_of_alignments().try_into().unwrap_or(255),
                            read.mapq(),
                        ),
                    };
                    res.push((AnnotatedRead::Counted(Box::new(info)), *org_index));
                    *org_index += 1;
                }
            }
        }
    }

    fn write_annotated_reads(
        bam: &mut rust_htslib::bam::IndexedReader,
        chunk: &Chunk,
        mut idx_to_annotated: HashMap<usize, AnnotatedRead>,
        out_bam_path: &Path,
        header: &rust_htslib::bam::Header,
        max_skip_len: u32,
        interner: &OurInterner,
    ) -> Result<()> {
        let mut out_bam = rust_htslib::bam::Writer::from_path(
            out_bam_path.join(format!("{}.bam", chunk.str_id())),
            header,
            rust_htslib::bam::Format::Bam,
        )?;
        /* let idx_to_annotated: HashMap<usize, &AnnotatedRead> = annotated_reads
        .iter()
        .map(|(read, org_index)| (*org_index, read))
        .collect(); */
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
            if let Some(anno_read) = idx_to_annotated.get_mut(&ii) {
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
                    AnnotatedRead::BarcodeNotInWhitelist(uncorrected_barcode) => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(4))?;
                        read.replace_aux(
                            b"CR",
                            rust_htslib::bam::record::Aux::String(
                                std::str::from_utf8(uncorrected_barcode).unwrap_or("non-utf8"),
                            ),
                        )?;
                    }
                    AnnotatedRead::NoBarcode => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(5))?;
                    }

                    AnnotatedRead::NoUMI => {
                        read.replace_aux(b"XF", rust_htslib::bam::record::Aux::U8(6))?;
                    }
                    AnnotatedRead::Counted(info) => {
                        //we have a read that was annotated
                        //write it to the output bam
                        let mut tag = String::new();
                        let mut first = true;

                        for gene in info.hits.correct.iter().sorted() {
                            if !first {
                                tag.push(',')
                            }
                            first = false;
                            tag.push_str(
                                interner.resolve(*gene).expect("string de-interning failed"),
                            );
                        }
                        read.replace_aux(b"XQ", rust_htslib::bam::record::Aux::String(&tag))?;

                        let mut tag = String::new();
                        let mut first = true;
                        for gene in info.hits.reverse.iter().sorted() {
                            if !first {
                                tag.push(',')
                            }
                            first = false;
                            tag.push_str(
                                interner.resolve(*gene).expect("string de-interning failed"),
                            );
                        }
                        read.replace_aux(b"XR", rust_htslib::bam::record::Aux::String(&tag))?;
                        read.replace_aux(
                            b"XP",
                            //convert back into sam's 1 based coordinates.
                            rust_htslib::bam::record::Aux::I32(info.corrected_position + 1i32),
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
                res.push(format!("{}.bam", chunk_name));
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
    fn generate_chunks(&self, bam: rust_htslib::bam::IndexedReader) -> Result<Vec<Chunk>> {
        let cg = ChunkedGenome::new(&self.reference_to_aggregation_trees, bam)?; // can't get the ParallelBridge to work with our lifetimes.
        Ok(cg.iter().collect())
    }

    fn hits(
        &self,
        chunk: &Chunk,
        read: &rust_htslib::bam::record::Record,
        interner: &mut OurInterner,
    ) -> Result<(
        Vec<string_interner::symbol::SymbolU32>,
        Vec<string_interner::symbol::SymbolU32>,
    )> {
        use crate::config::MatchDirection;
        let (tree, gene_ids) = self
            .reference_to_count_trees
            .get(&chunk.chr)
            .expect("Chr not found in trees");
        let blocks = read.blocks();
        if let crate::config::OverlapMode::Union = self.count_strategy.overlap {
            let mut gene_nos_seen_match = Vec::new();
            let mut gene_nos_seen_reverse = Vec::new();
            //todo: I don't like having this duplication
            for iv in blocks.iter() {
                if chunk.interval_outside(iv.0, iv.1) {
                    // if this block is outside of the region
                    // don't count it at all.
                    // if it is on a block boundary
                    // only count it for the left side.
                    // which is ok, since we place the blocks to the right
                    // of our intervals.
                    continue;
                }
                //todo: consider using either a bitset for the overlap range,
                //or no overlap range at all when doing Union.
                for r in tree.find(iv.0..iv.1) {
                    let entry = r.data();
                    let gene_no = entry.0;
                    let region_strand = entry.1;
                    let target = match (
                        &self.count_strategy.direction,
                        read.is_reverse(),
                        region_strand,
                    ) {
                        (MatchDirection::Forward, false, Strand::Forward) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Forward, false, Strand::Reverse) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Forward, true, Strand::Forward) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Forward, true, Strand::Reverse) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Forward, _, Strand::Unstranded) => {
                            &mut gene_nos_seen_match
                        }

                        (MatchDirection::Reverse, false, Strand::Forward) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Reverse, false, Strand::Reverse) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Reverse, true, Strand::Forward) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Reverse, true, Strand::Reverse) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Reverse, _, Strand::Unstranded) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Ignore, _, _) => &mut gene_nos_seen_match,
                    };
                    let gene_id = interner.get_or_intern(&gene_ids[gene_no as usize]);
                    if !target.iter().any(|x| *x == gene_id) {
                        // if we haven't seen this gene yet, add it
                        target.push(gene_id);
                    }
                }
            }
            for gg in [&mut gene_nos_seen_match, &mut gene_nos_seen_reverse] {
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
            Ok((gene_nos_seen_match, gene_nos_seen_reverse))
        } else {
            let mut gene_nos_seen_match =
                HashMap::<string_interner::symbol::SymbolU32, Vec<std::ops::Range<u32>>>::new();
            let mut gene_nos_seen_reverse =
                HashMap::<string_interner::symbol::SymbolU32, Vec<std::ops::Range<u32>>>::new();
            let mut bases_aligned = 0u32;
            for iv in blocks.iter() {
                bases_aligned += iv.1 - iv.0;
                if chunk.interval_outside(iv.0, iv.1) {
                    // if this block is outside of the region
                    // don't count it at all.
                    // if it is on a block boundary
                    // only count it for the left side.
                    // which is ok, since we place the blocks to the right
                    // of our intervals.
                    continue;
                }
                //todo: consider using either a bitset for the overlap range,
                //or no overlap range at all when doing Union.
                for r in tree.find(iv.0..iv.1) {
                    let entry = r.data();
                    let found_interval = r.interval();
                    let intersected_interval =
                        found_interval.start.max(iv.0)..found_interval.end.min(iv.1);

                    let gene_no = entry.0;
                    let region_strand = entry.1;
                    let target = match (
                        &self.count_strategy.direction,
                        read.is_reverse(),
                        region_strand,
                    ) {
                        (MatchDirection::Forward, false, Strand::Forward) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Forward, false, Strand::Reverse) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Forward, true, Strand::Forward) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Forward, true, Strand::Reverse) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Forward, _, Strand::Unstranded) => {
                            &mut gene_nos_seen_match
                        }

                        (MatchDirection::Reverse, false, Strand::Forward) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Reverse, false, Strand::Reverse) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Reverse, true, Strand::Forward) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Reverse, true, Strand::Reverse) => {
                            &mut gene_nos_seen_reverse
                        }
                        (MatchDirection::Reverse, _, Strand::Unstranded) => {
                            &mut gene_nos_seen_match
                        }
                        (MatchDirection::Ignore, _, _) => &mut gene_nos_seen_match,
                    };

                    let gene_id = interner.get_or_intern(&gene_ids[gene_no as usize]);
                    match target.entry(gene_id) {
                        std::collections::hash_map::Entry::Occupied(mut e) => {
                            e.get_mut().push(intersected_interval);
                        }
                        std::collections::hash_map::Entry::Vacant(e) => {
                            e.insert(vec![intersected_interval]);
                        }
                    }
                }
            }
            let gene_nos_seen_match = apply_count_strategy(
                &self.count_strategy,
                gene_nos_seen_match,
                bases_aligned as usize,
            );
            let gene_nos_seen_reverse = apply_count_strategy(
                &self.count_strategy,
                gene_nos_seen_reverse,
                bases_aligned as usize,
            );

            Ok((gene_nos_seen_match, gene_nos_seen_reverse))
        }
    }
}

pub struct TagMatcher {
    tag: [u8; 2],
}

impl TagMatcher {
    fn generate_chunks(&self, mut bam: rust_htslib::bam::IndexedReader) -> Result<Vec<Chunk>> {
        let keep_tids = chunked_genome::tids_with_reads(&mut bam)?;

        Ok(bam
            .header()
            .target_names()
            .iter()
            .enumerate()
            .filter_map(|(tid, name)| {
                let name = String::from_utf8(name.to_vec())
                    .expect("Failed to convert target name to string");
                let tid: u32 = tid.try_into().expect("tid should fit into u32");
                if keep_tids.iter().any(|x| *x == tid) {
                    Some(Chunk::new(
                        name,
                        tid,
                        0,
                        u32::try_from(bam.header().target_len(tid).unwrap_or(0))
                            .expect("reference length > u32 capacitiy."),
                    ))
                } else {
                    None
                }
            })
            .collect())
    }

    fn hits(
        &self,
        _chunk: &Chunk,
        read: &rust_htslib::bam::Record,
        interner: &mut OurInterner,
    ) -> Result<(
        Vec<string_interner::symbol::SymbolU32>,
        Vec<string_interner::symbol::SymbolU32>,
    )> {
        let mut genes_hit_correct = Vec::new();
        let genes_hit_reverse = Vec::new();
        if let Ok(rust_htslib::bam::record::Aux::String(value)) = read.aux(&self.tag) {
            genes_hit_correct.push(interner.get_or_intern(value));
        }
        Ok((genes_hit_correct, genes_hit_reverse))
    }
}

pub struct ReferenceMatcher {
    tids_to_include: Vec<i32>,
    direction: MatchDirection,
}

impl ReferenceMatcher {
    fn generate_chunks(&self, bam: rust_htslib::bam::IndexedReader) -> Result<Vec<Chunk>> {
        Ok(bam
            .header()
            .target_names()
            .iter()
            .enumerate()
            .filter_map(|(tid, name)| {
                let tid: u32 = tid.try_into().expect("tid should fit into u32");
                if self.tids_to_include.iter().any(|x| *x as u32 == tid) {
                    let name = std::str::from_utf8(*name)
                        .expect("Failed to convert target name to string")
                        .to_string();
                    Some(Chunk::new(
                        name,
                        tid as u32,
                        0,
                        u32::try_from(bam.header().target_len(tid).unwrap_or(0))
                            .expect("reference length > u32 capacitiy."),
                    ))
                } else {
                    None
                }
            })
            .collect())
    }

    fn hits(
        &self,
        chunk: &Chunk,
        read: &rust_htslib::bam::Record,
        interner: &mut OurInterner,
    ) -> Result<(
        Vec<string_interner::symbol::SymbolU32>,
        Vec<string_interner::symbol::SymbolU32>,
    )> {
        let mut genes_hit_correct = Vec::new();
        let mut genes_hit_reverse = Vec::new();
        if match (&self.direction, read.is_reverse()) {
            (MatchDirection::Ignore, _) => true,
            (MatchDirection::Forward, true) => false,
            (MatchDirection::Forward, false) => true,
            (MatchDirection::Reverse, true) => true,
            (MatchDirection::Reverse, false) => false,
        } {
            genes_hit_correct.push(interner.get_or_intern(&chunk.chr));
        }
        else 
        {
            genes_hit_reverse.push(interner.get_or_intern(&chunk.chr));
        }
        Ok((genes_hit_correct, genes_hit_reverse))
    }
}

fn merged_interval_length(ivs: &mut [std::ops::Range<u32>]) -> usize {
    if ivs.is_empty() {
        return 0;
    }
    ivs.sort_by(|a, b| a.start.cmp(&b.start).then_with(|| a.end.cmp(&b.end)));
    let mut current_start = ivs[0].start;
    let mut current_end = ivs[0].end;
    let mut total = 0;
    for iv in ivs.iter().skip(1) {
        if iv.start <= current_end {
            // overlap, merge
            current_end = current_end.max(iv.end);
        } else {
            // no overlap, push the current interval and start a new one
            total += (current_end - current_start) as usize;
            current_start = iv.start;
            current_end = iv.end;
        }
    }
    total += (current_end - current_start) as usize;
    total
}

fn apply_count_strategy(
    count_strategy: &crate::config::Strategy,
    gene_nos: HashMap<string_interner::symbol::SymbolU32, Vec<std::ops::Range<u32>>>,
    bases_aligned: usize,
) -> Vec<string_interner::symbol::SymbolU32> {
    //now merge the intervals, and convert them into length
    let mut gg_len = gene_nos
        .into_iter()
        .map(|(k, mut v)| (k, merged_interval_length(&mut v)))
        .collect::<HashMap<SymbolU32, usize>>();
    match count_strategy.overlap {
        crate::config::OverlapMode::Union => {
            unreachable!();
            // do nothing, we keep them all
        }
        crate::config::OverlapMode::IntersectionStrict => {
            //only keep those that are fully contained in the region
            gg_len.retain(|_, v| *v == bases_aligned);
        }
        crate::config::OverlapMode::IntersectionNonEmpty => {
            let any_fully_contained = gg_len.values().any(|v| *v == bases_aligned);
            if any_fully_contained {
                // only keep those that are fully contained in the region
                gg_len.retain(|_, v| *v == bases_aligned);
            } else {
                // multiple partial overlaps, keep all.
            }
        }
    }
    match count_strategy.multi_region {
        crate::config::MultiRegionHandling::Drop => {
            if gg_len.len() > 1 {
                // if there are multiple genes, drop them
                gg_len.clear();
            }
        }
        crate::config::MultiRegionHandling::CountBoth => {
            //do nothing.
        }
    }
    gg_len.into_keys().collect()
}
#[cfg(test)]
mod test {
    #[test]
    fn test_merged_interval_lengths() {
        use super::merged_interval_length;
        assert_eq!(merged_interval_length(&mut []), 0);
        assert_eq!(merged_interval_length(&mut [0..10, 10..20, 20..30]), 30);
        assert_eq!(merged_interval_length(&mut [0..10, 5..15, 10..20]), 20);
        assert_eq!(merged_interval_length(&mut [0..10, 5..15, 10..25]), 25);
        assert_eq!(merged_interval_length(&mut [0..10, 20..55]), 45);
        assert_eq!(
            merged_interval_length(&mut [20..30, 45..50, 28..32]),
            10 + 5 + 2
        );
    }
}
