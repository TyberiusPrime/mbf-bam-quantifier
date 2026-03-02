use anyhow::{Context, Result, bail};
use bstr::{BStr, BString};
use hamming_resonate::HammingResonatorWeighted;
use rust_htslib::htslib::htsExactFormat_unknown_format;
use std::{collections::{HashMap, HashSet}, path::PathBuf};

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
                Ok(HammingResonatorWeighted::with_max_dist(res, self.max_hamming as u32)
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
                    out.extend(best.as_bytes());
                }
                None => return None,
            }
        }
        Some(out)
    }

    pub fn correct_mtx_per_barcode(
        &self,
        feature_filename: &PathBuf,
        barcode_filename: &PathBuf,
        matrix_filename: &PathBuf,
    ) -> Result<()> {
        let raw_dir = matrix_filename
            .parent()
            .map(|p| p.join("raw"))
            .unwrap_or_else(|| PathBuf::from("."))
        ;
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

        //now the features stay the same.
        //so just symlink the file
        std::os::unix::fs::symlink(&raw_feature_filename, feature_filename).with_context(|| {
            format!(
                "Failed to symlink feature file: {:?} -> {:?}",
                raw_feature_filename, feature_filename
            )
        })?;

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
        let barcode_to_column: HashMap<&BStr, usize> = todo!();
        let matrix_coo = read_mtx(&raw_matrix_filename)?;
        //sum up the matrix in barcode direction.
        let barcode_counts = todo!();//column sums
        //now convert the existing HammingResonatorWeighted (which has all 0 scores)
        //into one where the scores == barcode counts.

        //then go through all the barcodes.
        //if the best matching in the HammingResonatorWeighted has a distance > 0,
        //we need to fold it's row into the one of the best match (add)
        //otherwise, we need to mark it as being kept.
        //
        //then filter it those being kept,
        //and write it out.
        //Also write out a new barcode file


        let mut keep = Vec::new();


        Ok(())
    }
}

fn read_mtx(filename: &PathBuf) -> Result<nalgebra_sparse::CooMatrix<i64>> {
    Ok(nalgebra_sparse::io::load_coo_from_matrix_market_file(filename)?)
}
