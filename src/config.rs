use anyhow::{Context, Result};
use noodles::bam;
use serde::{Deserialize, Serialize};

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
    pub fn get_bam_reader(
        &self,
    ) -> Result<bam::io::reader::Reader<noodles_bgzf::io::reader::Reader<std::fs::File>>> {
        Ok(bam::io::reader::Builder::default()
            .build_from_path(&self.bam)
            .with_context(|| format!("Failed to open bam file {} (without index)", &self.bam))?)
    }

    pub fn get_indexed_bam_reader(
        &self,
    ) -> Result<
        bam::io::indexed_reader::IndexedReader<noodles_bgzf::io::reader::Reader<std::fs::File>>,
    > {
        Ok(bam::io::indexed_reader::Builder::default()
            .build_from_path(&self.bam)
            .with_context(|| format!("Failed to open bam file {} (with index)", &self.bam))?)
    }
    /* fn get_bam_indexd_reader(&self) -> Result<u8> {
            Ok(bam::io::indexed_reader::Builder().default().build_from_path(&self.bam)?)
    } */

    pub fn get_bam_index(&self) -> Result<bam::bai::Index> {
        Ok(bam::bai::fs::read(format!("{}.bai", self.bam))
           .with_context(|| format!("Failed to open bai file {}.bai", &self.bam))?)
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
