use anyhow::{Context, Result};
use std::{collections::HashSet, path::PathBuf};

use crate::extractors::{self, UMIExtractor};
use serde::Deserializer;
use serde_valid::Validate;

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
fn validate_mode(mode: &extractors::UMIExtraction) -> Result<(), serde_valid::validation::Error> {
    if matches!(mode, extractors::UMIExtraction::NoUMI(_)) {
        Err(serde_valid::validation::Error::Custom(
            "Must define an extraction mode for barcodes".to_string(),
        ))
    } else {
        Ok(())
    }
}

type Whitelist = HashSet<Vec<u8>>;

#[derive(serde::Deserialize, Debug, Clone, Validate)]
#[serde(deny_unknown_fields)]
#[validate(custom = |s|validate_mode(&s.extract))]
pub struct CellBarcodes {
    extract: extractors::UMIExtraction,
    #[serde(deserialize_with = "u8_from_char_or_number")]
    separator_char: u8,
    #[serde(default)]
    max_hamming: u16,
    #[validate(min_length = 1)]
    whitelist_files: Vec<PathBuf>,

    #[serde(skip)]
    whitelists: Vec<Whitelist>,
}

impl CellBarcodes {
    pub fn init(&mut self) -> Result<()> {
        let wl: Result<_> = self
            .whitelist_files
            .iter()
            .map(|file| {
                Ok(std::fs::read_to_string(file)
                    .with_context(|| format!("Failed to read whitelist file: {:?}", file))?
                    .lines()
                    .map(|line| line.trim().as_bytes().to_vec())
                    .collect::<HashSet<_>>())
            })
            .collect();
        self.whitelists = wl?;
        Ok(())
    }

    pub fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<Vec<u8>> {
        self.extract.extract(read)
    }

    pub fn correct(&self, barcode: &[u8]) -> Option<Vec<u8>> {
        // possibly microopt: use cow...
        let parts = barcode.split(|&b| b == self.separator_char);
        let mut out = Vec::new();
        for (part, whitelist) in parts.zip(self.whitelists.iter()) {
            if whitelist.contains(part) {
                if !out.is_empty() {
                    out.push(self.separator_char);
                }
                out.extend(part);
            } else {
                match self.find_closest_by_hamming(part, whitelist) {
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
        whitelist: &'a Whitelist,
    ) -> Option<&'a [u8]> {
        use bio::alignment::distance::hamming;
        if self.max_hamming == 0 {
            return None; // No correction allowed
        }
        for entry in whitelist {
            if hamming(entry, part) <= self.max_hamming as u64 {
                return Some(entry);
            }
        }
        None
    }
}
