use anyhow::{Context, Result, bail};
use bstr::{BStr, BString, ByteSlice};
use hamming_resonate::HammingResonatorWeighted;
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    path::PathBuf,
};

use crate::extractors::{self, Extractors};
use serde::Deserializer;

pub fn u8_from_char_or_number<'de, D>(deserializer: D) -> Result<u8, D::Error>
where
    D: Deserializer<'de>,
{
    struct Visitor;

    impl serde::de::Visitor<'_> for Visitor {
        type Value = u8;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("either a byte character or a number 0..255")
        }

        fn visit_u64<E>(self, v: u64) -> Result<Self::Value, E>
        where
            E: serde::de::Error,
        {
            u8::try_from(v).map_err(|_| E::custom("Number too large for u8/char"))
        }

        fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
        where
            E: serde::de::Error,
        {
            match v.len() {
                0 => Err(E::custom("empty string")),
                1 => Ok(v.bytes().next().unwrap()),
                _ => Err(E::custom("string should be exactly one character long")),
            }
        }
    }

    deserializer.deserialize_any(Visitor)
}

type Whitelist = hamming_resonate::HammingResonatorWeighted;

#[derive(serde::Deserialize, Debug)]
pub(crate) enum AllowListMode {
    /// check each read against the allow list and correct if possible
    PerRead,
    /// First, count everything. Then correct 'virtual cells'
    /// to the closest barcode (within max_hamming distance), breaking ties towards
    /// barcodes with more counts.
    PerBarcode,
}

#[derive(serde::Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct CellBarcodes {
    extract: extractors::Extractor,
    #[serde(deserialize_with = "u8_from_char_or_number")]
    separator_char: u8, //todo: make optional?
    #[serde(default)]
    max_hamming: u16,
    #[serde(default)]
    allowlist_files: Vec<PathBuf>,
    pub(crate) allowlist_mode: AllowListMode,
    #[serde(skip)]
    pub(crate) allowlists: Vec<Whitelist>,
}

impl CellBarcodes {
    pub fn init(&mut self) -> Result<()> {
        let wl: Result<Vec<HammingResonatorWeighted>> = self
            .allowlist_files
            .iter()
            .map(|file| {
                let res = std::fs::read_to_string(file)
                    .with_context(|| format!("Failed to read allowlist file: {:?}", file))?
                    .lines()
                    .map(|line| (BString::from(line.trim().as_bytes()), 0.0f32))
                    .collect::<Vec<_>>();
                let lengths_observed: HashSet<usize> = res.iter().map(|x| x.0.len()).collect();
                if lengths_observed.len() >1 {
                    bail!("More than one length in allow list file. Barcodes must be of uniform length. Observed: {:?}", lengths_observed);
                }
                let deduped: HashSet<_> = res.iter().map(|x| &x.0).collect();
                if deduped.len() != res.len() {
                    bail!("Duplicate entries in allow list file: {:?}", file);
                }
                Ok(HammingResonatorWeighted::new(res, self.max_hamming as u32)
                    .with_context(|| format!("Failed to create HammingResonator for file: {:?}", file))?)
            })
            .collect();
        self.allowlists = wl?;
        Ok(())
    }

    pub fn check(&self, _config: &crate::config::Config) -> Result<()> {
        Ok(())
    }

    pub fn extract(&self, read: &rust_htslib::bam::record::Record) -> Result<Option<Vec<u8>>> {
        self.extract.extract(read)
    }

    /// Correct reads by hamming distance matched barcode,
    /// correcting each read individually while it's being processed.
    pub fn correct(&self, barcode: &[u8]) -> Option<Vec<u8>> {
        use bstr::ByteSlice;
        if barcode.is_empty() {
            return None;
        }
        if !matches!(self.allowlist_mode, AllowListMode::PerRead) {
            return Some(barcode.to_vec());
        }
        // possibly microopt: use cow...
        if self.allowlists.is_empty() {
            return Some(barcode.to_vec());
        }
        let parts = barcode.split(|&b| b == self.separator_char);
        let mut out = Vec::new();
        for (part, allowlist) in parts.zip(self.allowlists.iter()) {
            let best = allowlist.query_best(BStr::new(part)).ok()?; // length mismatch -> no match
            match best {
                Some(best) => {
                    if !out.is_empty() {
                        out.push(self.separator_char);
                    }
                    out.extend(best.0.as_bytes());
                }
                None => return None,
            }
        }
        //since we queried to be within hamming distance on each part
        //we still need to verify that the total is inside the max hamming distance
        if self.allowlists.len() > 1
            && hamming_resonate::hamming_distance(barcode, &out) > self.max_hamming as u32
        {
            return None;
        }
        Some(out)
    }

    fn rename_raw_files(
        &self,
        feature_filename: &PathBuf,
        barcode_filename: &PathBuf,
        matrix_filename: &PathBuf,
        matrix_stats_filename: &PathBuf,
    ) -> Result<(PathBuf, PathBuf)> {
        let raw_dir = matrix_filename
            .parent()
            .map(|p| p.join("raw"))
            .unwrap_or_else(|| PathBuf::from("."));
        std::fs::create_dir_all(&raw_dir)
            .with_context(|| format!("Failed to create raw directory: {:?}", raw_dir))?;
        //now let's move all three files into the raw dir
        let raw_feature_filename = raw_dir.join(feature_filename.file_name().unwrap());
        std::fs::rename(feature_filename, &raw_feature_filename).with_context(|| {
            format!(
                "Failed to move feature file to raw directory: {:?} -> {:?}",
                feature_filename, raw_feature_filename
            )
        })?;
        let raw_barcode_filename = raw_dir.join(barcode_filename.file_name().unwrap());
        std::fs::rename(barcode_filename, &raw_barcode_filename).with_context(|| {
            format!(
                "Failed to move barcode file to raw directory: {:?} -> {:?}",
                barcode_filename, raw_barcode_filename
            )
        })?;
        let raw_matrix_filename = raw_dir.join(matrix_filename.file_name().unwrap());
        std::fs::rename(matrix_filename, &raw_matrix_filename).with_context(|| {
            format!(
                "Failed to move matrix file to raw directory: {:?} -> {:?}",
                matrix_filename, raw_matrix_filename
            )
        })?;

        let raw_matrix_stats_filename = raw_dir.join(matrix_stats_filename.file_name().unwrap());
        std::fs::rename(matrix_stats_filename, &raw_matrix_stats_filename).with_context(|| {
            format!(
                "Failed to move matrix stats file to raw directory: {:?} -> {:?}",
                matrix_stats_filename, raw_matrix_stats_filename
            )
        })?;

        //now the features stay the same.
        //so just symlink the file
        std::os::unix::fs::symlink(&raw_feature_filename, feature_filename).with_context(|| {
            format!(
                "Failed to symlink feature file: {:?} -> {:?}",
                raw_feature_filename, feature_filename
            )
        })?;
        Ok((raw_barcode_filename, raw_matrix_filename))
    }

    fn read_matrix_and_quantify_barcodes(
        &self,
        raw_barcode_filename: PathBuf,
        raw_matrix_filename: PathBuf,
    ) -> Result<(
        Vec<BString>,
        HammingResonatorWeighted,
        nalgebra_sparse::CooMatrix<i64>,
        Vec<i64>,

    )> {
        use itertools::Itertools;
        //read the barcodes, get us a hamming resonator with
        //scores =  counts
        let barcodes: Vec<BString> = {
            //read all the lines
            std::fs::read_to_string(&raw_barcode_filename)
                .with_context(|| {
                    format!("Failed to read barcode file: {:?}", raw_barcode_filename)
                })?
                .lines()
                .map(|line| BString::from(line.as_bytes()))
                .collect()
        };
        let barcode_to_column: HashMap<&BStr, usize> = barcodes
            .iter()
            .enumerate()
            .map(|(i, b)| (b.as_ref(), i))
            .collect();
        let matrix_coo =
            nalgebra_sparse::io::load_coo_from_matrix_market_file(&raw_matrix_filename)?;
        assert_eq!(
            barcodes.len(),
            matrix_coo.ncols(),
            "Number of barcodes must match number of columns in matrix"
        );
        // Column sums: total counts per barcode column.
        let mut barcode_sums = vec![0i64; barcodes.len()];
        for (&col, &val) in matrix_coo
            .col_indices()
            .iter()
            .zip(matrix_coo.values().iter())
        {
            barcode_sums[col] += val;
        }

        //this is still pre-combinatorial.
        let allowed_barcodes: Vec<Vec<BString>> =
            self.allowlists.iter().map(|x| x.to_seqs()).collect();
        let mut allowed_barcodes_scored: Vec<(BString, f32)> = allowed_barcodes
            .iter()
            .multi_cartesian_product()
            .map(|parts| {
                BString::from(
                    parts
                        .iter()
                        .map(|s| s.as_bytes())
                        .collect::<Vec<&[u8]>>()
                        .join(&self.separator_char),
                )
            })
            .map(|barcode| {
                let score = match barcode_to_column.get(barcode.as_bstr()) {
                    Some(&col) => barcode_sums[col] as f32,
                    None => 0.0,
                };
                (barcode, score)
            })
            .collect();
        allowed_barcodes_scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap()); //sort
        //alphabetically
        let hamming_db =
            HammingResonatorWeighted::new(allowed_barcodes_scored, self.max_hamming as u32)
                .context("Failed to build HammingResonator for corrected barcodes")?;

        Ok((barcodes, hamming_db, matrix_coo, barcode_sums))
    }

    fn write_corrected_matrix(
        new_coo: nalgebra_sparse::CooMatrix<i64>,
        hamming_db: &HammingResonatorWeighted,
        barcode_filename: &PathBuf,
        matrix_filename: &PathBuf,
    ) -> Result<()> {
        // Converting to CSC automatically sums duplicate (row, col) entries.
        let csc = nalgebra_sparse::CscMatrix::from(&new_coo);

        // Write the corrected matrix in MatrixMarket format.
        {
            nalgebra_sparse::io::save_to_matrix_market_file(&csc, matrix_filename)?;
        }

        // Write the complete barcode file.
        {
            use std::io::Write;
            let mut bw = std::io::BufWriter::new(
                std::fs::File::create(barcode_filename).with_context(|| {
                    format!(
                        "Failed to create corrected barcodes: {:?}",
                        barcode_filename
                    )
                })?,
            );
            for bc in hamming_db.to_seqs() {
                bw.write_all(bc.as_slice())?;
                bw.write_all(b"\n")?;
            }
        }
        Ok(())
    }

    fn write_correction_stats(
        matrix_stats_filename: &PathBuf,
        total_perfect: i64,
        total_corrected: i64,
        total_uncorrectable: i64,
    ) -> Result<()> {
        use std::io::Write;
        let mut bw = std::io::BufWriter::new(
            std::fs::File::create(matrix_stats_filename).with_context(|| {
                format!(
                    "Failed to create matrix stats file: {:?}",
                    matrix_stats_filename
                )
            })?,
        );
        writeln!(bw, "stat\tcount\n")?;
        writeln!(bw, "total_perfect\t{}", total_perfect)?;
        writeln!(bw, "total_corrected\t{}", total_corrected)?;
        writeln!(bw, "total_uncorrectable\t{}", total_uncorrectable)?;
        Ok(())
    }

    /// Correct barcodes after counting everything.
    /// We match them to the nearest by hamming distance,
    /// and fold them onto the 'perfect' barcodes,
    /// breaking ties by choosing the higher counted perfect barcode
    pub fn correct_mtx_per_barcode(
        &self,
        feature_filename: &PathBuf,
        barcode_filename: &PathBuf,
        matrix_filename: &PathBuf,
        matrix_stats_filename: &PathBuf,
    ) -> Result<()> {
        assert!(matches!(self.allowlist_mode, AllowListMode::PerBarcode));
        assert!(!self.allowlist_files.is_empty());
        let (raw_barcode_filename, raw_matrix_filename) = self.rename_raw_files(
            feature_filename,
            barcode_filename,
            matrix_filename,
            matrix_stats_filename,
        )?;

        let (observed_barcodes_in_matrix_order, hamming_db, matrix_coo, barcode_sums) =
            self.read_matrix_and_quantify_barcodes(raw_barcode_filename, raw_matrix_filename)?;

        let barcode_to_column: BTreeMap<BString, usize> = hamming_db
            .to_seqs()
            .into_iter()
            .enumerate()
            .map(|(i, b)| (b, i))
            .collect();

        let mut col_remap: Vec<Option<usize>> = vec![];
        //we are going to make use of the fact that the CSC
        //constructor sums up triplets
        let mut total_perfect = 0;
        let mut total_corrected = 0;
        let mut total_uncorrectable = 0;
        for (barcode_col, barcode) in observed_barcodes_in_matrix_order.iter().enumerate() {
            if let Some((best_barcode, hamming_dist, _count)) =
                hamming_db.query_best(barcode.as_bstr())?
            {
                col_remap.push(Some(*barcode_to_column.get(best_barcode).unwrap()));
                if hamming_dist == 0 {
                    total_perfect += barcode_sums[barcode_col];
                }else {
                    total_corrected += barcode_sums[barcode_col];
                }
            } else {
                col_remap.push(None);
                total_uncorrectable += barcode_sums[barcode_col];
            }
        }

        // Build a new COO matrix with remapped columns (duplicates arise when multiple
        // raw barcodes fold into the same target).
        let n_features = matrix_coo.nrows();
        let n_new_cols = barcode_to_column.len();
        let mut new_coo = nalgebra_sparse::CooMatrix::<i64>::new(n_features, n_new_cols);
        for ((&row, &old_col), &val) in matrix_coo
            .row_indices()
            .iter()
            .zip(matrix_coo.col_indices().iter())
            .zip(matrix_coo.values().iter())
        {
            if let Some(new_col) = col_remap[old_col] {
                new_coo.push(row, new_col, val);
            }
        }
        Self::write_corrected_matrix(new_coo, &hamming_db, barcode_filename, matrix_filename)?;
        Self::write_correction_stats(matrix_stats_filename, total_perfect, total_corrected, total_uncorrectable)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::path::Path;

    fn write_file(path: &Path, content: &str) {
        std::fs::write(path, content).unwrap();
    }

    fn read_lines(path: &Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.to_string())
            .collect()
    }

    fn parse_stats_file(path: &Path) -> HashMap<String, i64> {
        let content = std::fs::read_to_string(path).unwrap();
        let mut stats = HashMap::new();
        for line in content.lines().skip(1) { // Skip header line
            if line.trim().is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() == 2 {
                let stat_name = parts[0].to_string();
                let stat_value = parts[1].parse::<i64>().unwrap();
                stats.insert(stat_name, stat_value);
            }
        }
        stats
    }

    /// Parse a corrected MTX + barcodes file into a map of (barcode_name, row_0based) -> count.
    fn parse_corrected(matrix_path: &Path, barcode_path: &Path) -> HashMap<(String, usize), i64> {
        let barcodes = read_lines(barcode_path);
        let matrix_coo =
            nalgebra_sparse::io::load_coo_from_matrix_market_file(matrix_path).unwrap();
        let mut result = HashMap::new();
        for (i, j, v) in matrix_coo.triplet_iter() {
            let bc = barcodes[j].clone();
            let key = (bc, i);
            match result.entry(key) {
                std::collections::hash_map::Entry::Occupied(_) => panic!("duplicate MTX entry"),
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(*v);
                }
            }
        }
        for gene_idx in 0..matrix_coo.nrows() {
            for bc in &barcodes {
                result.entry((bc.clone(), gene_idx)).or_insert(0);
            }
        }
        result
    }

    /// Run PerBarcode correction in a fresh temp directory.
    ///
    /// Returns `(dir, sorted_barcodes, counts_map, stats)`. The caller must hold `dir`
    /// for as long as file-system checks are needed; it is cleaned up on drop.
    fn run_per_barcode_correction(
        allowlist_content: &str,
        features_content: &str,
        barcodes_content: &str,
        matrix_content: &str,
        max_hamming: u16,
    ) -> (
        tempfile::TempDir,
        Vec<String>,
        HashMap<(String, usize), i64>,
        HashMap<String, i64>,
    ) {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path();

        let allowlist = p.join("allowlist.txt");
        let features = p.join("features.tsv");
        let barcodes_file = p.join("barcodes.tsv");
        let matrix_file = p.join("matrix.mtx");
        let matrix_stats_file = p.join("matrix.mtx.stats.tsv");

        write_file(&allowlist, allowlist_content);
        write_file(&features, features_content);
        write_file(&barcodes_file, barcodes_content);
        write_file(&matrix_file, matrix_content);
        write_file(&matrix_stats_file, "left blank"); // content not relevant for this test

        let mut cb = CellBarcodes {
            extract: crate::extractors::Extractor::Tag(crate::extractors::Tag {
                tag: [b'C', b'B'],
            }),
            separator_char: b'|',
            max_hamming,
            allowlist_files: vec![allowlist],
            allowlist_mode: AllowListMode::PerBarcode,
            allowlists: Vec::new(),
        };
        cb.init().unwrap();
        cb.correct_mtx_per_barcode(&features, &barcodes_file, &matrix_file, &matrix_stats_file)
            .unwrap();

        let counts = parse_corrected(&matrix_file, &barcodes_file);
        let barcodes = read_lines(&barcodes_file);
        let stats = parse_stats_file(&matrix_stats_file);
        (dir, barcodes, counts, stats)
    }

    // ── basic ──────────────────────────────────────────────────────────────────

    /// Two valid 4-mers (AAAA, CCCC) plus one unused allowlist entry (TTTT).
    /// One noisy variant of each valid barcode is present in the observed data;
    /// a completely off-target barcode (GGGG, hamming 4 from everything) is dropped.
    ///
    /// Barcodes and expected fold-in:
    ///   AAAA (valid)   total counts 6+4 = 10
    ///   AAAT (noisy)   total counts 2+1 =  3  → folds into AAAA
    ///   CCCC (valid)   total counts 15+5 = 20
    ///   CCCG (noisy)   total counts 3        → folds into CCCC
    ///   GGGG (off-target, hamming ≥4)        → dropped
    #[test]
    fn test_correct_mtx_per_barcode_basic() {
        let (dir, barcodes, counts, stats) = run_per_barcode_correction(
            "TTTT\nAAAA\nCCCC\n",
            "geneA\ngeneB\n",
            "AAAA\nAAAT\nCCCC\nCCCG\nGGGG\n",
            "%%MatrixMarket matrix coordinate integer general\n\
             %\n\
             2 5 9\n\
             1 1 6\n\
             2 1 4\n\
             1 2 2\n\
             2 2 1\n\
             1 3 15\n\
             2 3 5\n\
             1 4 3\n\
             1 5 30\n\
             2 5 3\n",
            1,
        );

        // raw/ must exist and contain all three original files.
        let raw = dir.path().join("raw");
        assert!(
            raw.join("features.tsv").exists(),
            "raw/features.tsv missing"
        );
        assert!(
            raw.join("barcodes.tsv").exists(),
            "raw/barcodes.tsv missing"
        );
        assert!(raw.join("matrix.mtx").exists(), "raw/matrix.mtx missing");
        assert!(
            raw.join("matrix.mtx.stats.tsv").exists(),
            "raw/matrix.mtx.stats.tsv missing"
        );

        // features.tsv must be a symlink pointing into raw/.
        assert!(
            dir.path()
                .join("features.tsv")
                .symlink_metadata()
                .unwrap()
                .file_type()
                .is_symlink(),
            "features.tsv should be a symlink"
        );

        // Output must contain exactly the three allowlist entries.
        assert_eq!(barcodes, vec!["AAAA", "CCCC", "TTTT"]);

        // Folded counts.
        assert_eq!(counts[&("AAAA".into(), 0)], 8, "AAAA geneA = 6+2");
        assert_eq!(counts[&("AAAA".into(), 1)], 5, "AAAA geneB = 4+1");
        assert_eq!(counts[&("CCCC".into(), 0)], 18, "CCCC geneA = 15+3");
        assert_eq!(counts[&("CCCC".into(), 1)], 5, "CCCC geneB");
        assert_eq!(counts[&("TTTT".into(), 0)], 0, "TTTT geneA");
        assert_eq!(counts[&("TTTT".into(), 1)], 0, "TTTT geneB");

        // Verify correction stats.
        // AAAA: perfect (10 total), CCCC: perfect (20 total) = 30 perfect
        // AAAT: corrected to AAAA (3 total), CCCG: corrected to CCCC (3 total) = 6 corrected  
        // GGGG: uncorrectable (33 total) = 33 uncorrectable
        assert_eq!(stats[&"total_perfect".to_string()], 30, "total_perfect");
        assert_eq!(stats[&"total_corrected".to_string()], 6, "total_corrected");
        assert_eq!(stats[&"total_uncorrectable".to_string()], 33, "total_uncorrectable");
    }

    // ── allowlist order invariance ─────────────────────────────────────────────

    /// The order of entries inside the allowlist file must not affect which
    /// barcode a noisy read is corrected to, nor the final counts.
    /// (No ties exist in this dataset, so tie-breaking index is irrelevant.)
    #[test]
    fn test_correct_mtx_per_barcode_allowlist_order_invariant() {
        let matrix = "%%MatrixMarket matrix coordinate integer general\n\
                      %\n\
                      2 5 9\n\
                      1 1 6\n\
                      2 1 4\n\
                      1 2 2\n\
                      2 2 1\n\
                      1 3 15\n\
                      2 3 5\n\
                      1 4 3\n\
                      1 5 30\n\
                      2 5 3\n";
        let barcodes_in = "AAAA\nAAAT\nCCCC\nCCCG\nGGGG\n";
        let features = "geneA\ngeneB\n";

        let run = |allowlist: &str| {
            let (_dir, barcodes, counts, stats) =
                run_per_barcode_correction(allowlist, features, barcodes_in, matrix, 1);
            (barcodes, counts, stats)
        };

        let (b1, c1, s1) = run("TTTT\nAAAA\nCCCC\n");
        let (b2, c2, s2) = run("CCCC\nTTTT\nAAAA\n");
        let (b3, c3, s3) = run("AAAA\nCCCC\nTTTT\n");

        assert_eq!(b1, b2, "barcodes differ for order 2");
        assert_eq!(b1, b3, "barcodes differ for order 3");
        assert_eq!(c1, c2, "counts differ for order 2");
        assert_eq!(c1, c3, "counts differ for order 3");
        assert_eq!(s1, s2, "stats differ for order 2");
        assert_eq!(s1, s3, "stats differ for order 3");
    }

    // ── larger hamming distance ────────────────────────────────────────────────

    /// With max_hamming=2:
    ///   AATT (hamming 2 from AAAA, 4 from CCCC) → folds into AAAA
    ///   GCGT (hamming 3 from CCCC, 4 from AAAA) → outside budget, dropped
    #[test]
    fn test_correct_mtx_per_barcode_hamming2() {
        let (_, barcodes, counts, stats) = run_per_barcode_correction(
            "AAAA\nCCCC\n",
            "geneA\n",
            "AAAA\nAATT\nCCCC\nGCGT\n",
            "%%MatrixMarket matrix coordinate integer general\n\
             %\n\
             1 4 4\n\
             1 1 10\n\
             1 2 3\n\
             1 3 20\n\
             1 4 5\n",
            2,
        );

        assert_eq!(barcodes, vec!["AAAA", "CCCC"]);
        // AATT (hamming 2) folds into AAAA.
        assert_eq!(counts[&("AAAA".into(), 0)], 13, "AAAA geneA = 10+3");
        // GCGT (hamming 3) is dropped entirely.
        assert_eq!(counts[&("CCCC".into(), 0)], 20, "CCCC geneA");

        // Verify correction stats.
        // AAAA: perfect (10), CCCC: perfect (20) = 30 perfect
        // AATT: corrected to AAAA (3) = 3 corrected
        // GCGT: uncorrectable (9 total = 4+5) = 9 uncorrectable
        assert_eq!(stats[&"total_perfect".to_string()], 30, "total_perfect");
        assert_eq!(stats[&"total_corrected".to_string()], 3, "total_corrected");
        assert_eq!(stats[&"total_uncorrectable".to_string()], 5, "total_uncorrectable");
    }

    // ── count-based tie-breaking ───────────────────────────────────────────────

    /// When a noisy barcode is equidistant (hamming 1) from two valid barcodes,
    /// it must fold into the one whose observed column has the higher total count.
    ///
    ///   AACC  (valid, col 1): geneA = 100
    ///   TTCC  (valid, col 2): geneA =  10
    ///   ATCC  (noisy): hamming 1 from AACC *and* TTCC → must fold into AACC (higher count)
    #[test]
    fn test_correct_mtx_per_barcode_count_preference() {
        let (_, barcodes, counts, stats) = run_per_barcode_correction(
            "AACC\nTTCC\n",
            "geneA\n",
            "AACC\nTTCC\nATCC\n",
            "%%MatrixMarket matrix coordinate integer general\n\
             %\n\
             1 3 3\n\
             1 1 100\n\
             1 2 10\n\
             1 3 5\n",
            1,
        );

        assert_eq!(barcodes, vec!["AACC", "TTCC"]);
        // ATCC folds into AACC (higher count wins the tie).
        assert_eq!(counts[&("AACC".into(), 0)], 105, "AACC geneA = 100+5");
        assert_eq!(counts[&("TTCC".into(), 0)], 10, "TTCC geneA");

        // Verify correction stats.
        // AACC: perfect (100), TTCC: perfect (10) = 110 perfect
        // ATCC: corrected to AACC (8 total = 3+5) = 8 corrected
        // No uncorrectable barcodes = 0 uncorrectable
        assert_eq!(stats[&"total_perfect".to_string()], 110, "total_perfect");
        assert_eq!(stats[&"total_corrected".to_string()], 5, "total_corrected");
        assert_eq!(stats[&"total_uncorrectable".to_string()], 0, "total_uncorrectable");
    }
    #[test]
    fn test_correct_mtx_per_barcode_count_preference_swap() {
        let (_, barcodes, counts, stats) = run_per_barcode_correction(
            "TTCC\nAACC\n",
            "geneA\n",
            "AACC\nTTCC\nATCC\n",
            "%%MatrixMarket matrix coordinate integer general\n\
             %\n\
             1 3 3\n\
             1 1 100\n\
             1 2 10\n\
             1 3 5\n",
            1,
        );

        assert_eq!(barcodes, vec!["AACC", "TTCC"]);
        // ATCC folds into AACC (higher count wins the tie).
        assert_eq!(counts[&("AACC".into(), 0)], 105, "AACC geneA = 100+5");
        assert_eq!(counts[&("TTCC".into(), 0)], 10, "TTCC geneA");

        // Verify correction stats (same as previous test - allowlist order shouldn't matter).
        // AACC: perfect (100), TTCC: perfect (10) = 110 perfect
        // ATCC: corrected to AACC (8 total = 3+5) = 8 corrected
        // No uncorrectable barcodes = 0 uncorrectable
        assert_eq!(stats[&"total_perfect".to_string()], 110, "total_perfect");
        assert_eq!(stats[&"total_corrected".to_string()], 5, "total_corrected");
        assert_eq!(stats[&"total_uncorrectable".to_string()], 0, "total_uncorrectable");
    }
}
