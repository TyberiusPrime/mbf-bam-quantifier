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
    )> {
        use itertools::Itertools;
        //read the barcodes, get us a hammng resonator with
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
        let allowed_barcodes_scored: Vec<(BString, f32)> = allowed_barcodes
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
        let hamming_db =
            HammingResonatorWeighted::new(allowed_barcodes_scored, self.max_hamming as u32)
                .context("Failed to build HammingResonator for corrected barcodes")?;

        Ok((barcodes, hamming_db, matrix_coo))
    }

    pub fn correct_mtx_per_barcode(
        &self,
        feature_filename: &PathBuf,
        barcode_filename: &PathBuf,
        matrix_filename: &PathBuf,
    ) -> Result<()> {
        assert!(matches!(self.allowlist_mode, AllowListMode::PerBarcode));
        assert!(!self.allowlist_files.is_empty());
        let (raw_barcode_filename, raw_matrix_filename) =
            self.rename_raw_files(feature_filename, barcode_filename, matrix_filename)?;

        let (observed_barcodes_in_matrix_order, hamming_db, matrix_coo) =
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
        for barcode in observed_barcodes_in_matrix_order.iter() {
            if let Some((best_barcode, _hamming_dist, _count)) =
                hamming_db.query_best(barcode.as_bstr())?
            {
                col_remap.push(Some(*barcode_to_column.get(best_barcode).unwrap()));
            } else {
                col_remap.push(None);
            }
        }
        dbg!(&observed_barcodes_in_matrix_order);
        dbg!(&col_remap);

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
                std::collections::hash_map::Entry::Occupied(_occupied_entry) => panic!(),
                std::collections::hash_map::Entry::Vacant(vacant_entry) => {
                    vacant_entry.insert(*v);
                }
            }
        }
        for gene_idx in 0..matrix_coo.nrows() {
            for bc in &barcodes {
                let key = (bc.clone(), gene_idx);
                result.entry(key).or_insert(0);
            }
        }
        result
    }

    #[test]
    fn test_correct_mtx_per_barcode_folds_noisy_barcodes() {
        let dir = tempfile::TempDir::new().unwrap();
        let p = dir.path();

        // Allowlist: two valid 4-mers.
        let allowlist = p.join("allowlist.txt");
        write_file(&allowlist, "TTTT\nAAAA\nCCCC\n");

        // features.tsv: two genes (row 1 = geneA, row 2 = geneB in MTX 1-based indexing).
        let features = p.join("features.tsv");
        write_file(&features, "geneA\ngeneB\n");

        // barcodes.tsv: two valid barcodes and one noisy variant each.
        //   col 1 = AAAA (valid),  total counts = 6+4 = 10
        //   col 2 = AAAT (noisy),  total counts = 2+1 =  3  → should fold into AAAA
        //   col 3 = CCCC (valid),  total counts = 15+5 = 20
        //   col 4 = CCCG (noisy),  total counts = 3   =  3  → should fold into CCCC
        let barcodes_file = p.join("barcodes.tsv");
        write_file(&barcodes_file, "AAAA\nAAAT\nCCCC\nCCCG\nGGGG\n");

        // matrix.mtx: 2 features × 4 barcodes, 7 non-zero entries.
        let matrix_file = p.join("matrix.mtx");
        write_file(
            &matrix_file,
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
        );

        let mut cb = CellBarcodes {
            extract: crate::extractors::Extractor::Tag(crate::extractors::Tag {
                tag: [b'C', b'B'],
            }),
            separator_char: b'|',
            max_hamming: 1,
            allowlist_files: vec![allowlist],
            allowlist_mode: AllowListMode::PerBarcode,
            allowlists: Vec::new(),
        };
        cb.init().unwrap();

        cb.correct_mtx_per_barcode(&features, &barcodes_file, &matrix_file)
            .unwrap();

        // raw/ must exist and contain all three original files.
        let raw = p.join("raw");
        assert!(
            raw.join("features.tsv").exists(),
            "raw/features.tsv missing"
        );
        assert!(
            raw.join("barcodes.tsv").exists(),
            "raw/barcodes.tsv missing"
        );
        assert!(raw.join("matrix.mtx").exists(), "raw/matrix.mtx missing");

        // features.tsv must be a symlink pointing into raw/.
        assert!(
            features
                .symlink_metadata()
                .unwrap()
                .file_type()
                .is_symlink(),
            "features.tsv should be a symlink"
        );

        // Corrected barcodes file must contain exactly the two valid ones.
        let out_barcodes = read_lines(&barcodes_file);
        assert_eq!(out_barcodes.len(), 3); //includes T
        assert!(out_barcodes.contains(&"AAAA".to_string()));
        assert!(out_barcodes.contains(&"CCCC".to_string()));
        assert!(out_barcodes.contains(&"TTTT".to_string()));

        // Corrected matrix counts:
        //   AAAA: geneA = 6+2 = 8,  geneB = 4+1 = 5
        //   CCCC: geneA = 15+3 = 18, geneB = 5
        //we have rewritten the matrix file.
        let bc = std::fs::read_to_string(&barcodes_file).unwrap();
        dbg!(bc);
        let counts = parse_corrected(&matrix_file, &barcodes_file);
        dbg!(&counts);
        assert_eq!(counts[&("AAAA".into(), 0)], 2 + 6, "AAAA geneA");
        assert_eq!(counts[&("AAAA".into(), 1)], 5, "AAAA geneB");
        assert_eq!(counts[&("CCCC".into(), 0)], 18, "CCCC geneA");
        assert_eq!(counts[&("CCCC".into(), 1)], 5, "CCCC geneB");

        assert_eq!(counts[&("TTTT".into(), 0)], 0, "TTTT geneA");
        assert_eq!(counts[&("TTTT".into(), 1)], 0, "TTTT geneB");
    }
}
