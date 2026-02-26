use anyhow::{Context, Result};
use anyhow::{bail, Context, Result};
use std::{collections::HashSet, path::PathBuf};

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

type Whitelist = HashSet<Vec<u8>>;

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct CellBarcodes {
    extract: extractors::Extractor,
    #[serde(deserialize_with = "u8_from_char_or_number")]
    separator_char: u8,  //todo: make optional?
    separator_char: u8, //todo: make optional?
    #[serde(default)]
    max_hamming: u16,
    #[serde(default)]
    allowlist_files: Vec<PathBuf>,

    #[serde(skip)]
    allowlists: Vec<Whitelist>,
}

impl CellBarcodes {
    pub fn init(&mut self) -> Result<()> {
        let wl: Result<_> = self
            .allowlist_files
            .iter()
            .map(|file| {
                let res = std::fs::read_to_string(file)
                    .with_context(|| format!("Failed to read allowlist file: {:?}", file))?
                    .lines()
                    .map(|line| line.trim().as_bytes().to_vec())
                    .collect::<HashSet<_>>();
                let lengths_observed: HashSet<usize> = res.iter().map(|x| x.len()).collect();
                if lengths_observed.len() >1 {
                    bail!("More than one length in allow list file. Barcodes must be of uniform length. Observed: {:?}", lengths_observed);
                }
                Ok(res)
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
        // possibly microopt: use cow...
        if self.allowlists.is_empty() {
            if barcode.is_empty() {
                return None;
            } else {
                return Some(barcode.to_vec());
            }
        }
        let parts = barcode.split(|&b| b == self.separator_char);
        let mut out = Vec::new();
        for (part, allowlist) in parts.zip(self.allowlists.iter()) {
            if allowlist.contains(part) {
                if !out.is_empty() {
                    out.push(self.separator_char);
                }
                out.extend(part);
            } else {
                match self.find_closest_by_hamming(part, allowlist) {
                    Some(corrected) => {
                        if !out.is_empty() {
                            out.push(self.separator_char);
                        }
                        out.extend(corrected);
                    }
                    None => return None,
                }
            }
        }
        Some(out)
    }

    fn find_closest_by_hamming<'a>(
        &self,
        part: &[u8],
        allowlist: &'a Whitelist,
    ) -> Option<&'a [u8]> {
        use bio::alignment::distance::hamming;
        if self.max_hamming == 0 {
            return None; // No correction allowed
        }
        allowlist
            .iter()
            .find(|&entry| hamming(entry, part) <= self.max_hamming as u64)
            .map(|v| v.as_slice())
    }
}
