use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::io;
use crate::quantification::{Quant, Quantification};

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Config {
    pub input: Input,
    #[serde(alias = "quant")]
    pub quantification: Quantification,
    pub output: Output,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]

pub struct Input {
    pub bam: String,
    pub gtf: Option<String>,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]

pub struct Output {
    pub directory: String,
}

impl Config {
    pub fn check(&self) -> Result<()> {
        self.quantification.check(&self)?;
        Ok(())
    }
}

impl Input {
    pub fn get_bam_reader(&self) -> Result<rust_htslib::bam::Reader> {
        Ok(rust_htslib::bam::Reader::from_path(&self.bam)
            .with_context(|| format!("Failed to open bam file {} (without index)", &self.bam))?)
    }

    pub fn get_indexed_bam_reader(&self) -> Result<rust_htslib::bam::IndexedReader> {
        Ok(rust_htslib::bam::IndexedReader::from_path(&self.bam)
            .with_context(|| format!("Failed to open bam file {} (with index)", &self.bam))?)
    }

    pub fn read_gtf(
        &self,
        accepted_features: &[&str],
    ) -> Result<HashMap<String, crate::gtf::GTFEntrys>> {
        crate::gtf::parse_ensembl_gtf(
            self.gtf
                .as_ref()
                .context("No GTF defined in input, but required")?,
            accepted_features
                .iter()
                .map(|s| s.to_string())
                .collect::<HashSet<_>>(),
        )
    }
    /* reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        if ignore_unmapped && record.reference_sequence_id().is_none() {
            continue;
        }

        if let Some(name) = record.name() {
            func(name);
        }
    } */
}
