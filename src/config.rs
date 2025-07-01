use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
};

use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

use crate::quantification::{Quant, Quantification};

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Config {
    pub input: Input,
    #[serde(default)]
    pub filter: Vec<crate::filters::Filter>,
    #[serde(alias = "quant")]
    pub quantification: Quantification,
    pub output: Output,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub bam: String,
    pub source: Source,
}

fn default_aggr_feature() -> String {
    "gene".to_string()
}

fn default_aggr_id_attribute() -> String {
    "gene_id".to_string()
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
#[serde(tag = "type")]
pub enum Source {
    #[serde(alias = "gtf")]
    GTF(GTFConfig),
    #[serde(alias = "bam_references")]
    BamReferences,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GTFConfig {
    pub filename: String,
    pub feature: String,
    pub id_attribute: String,
    #[serde(default = "default_aggr_feature")]
    pub aggr_feature: String,
    #[serde(default = "default_aggr_id_attribute")]
    pub aggr_id_attribute: String,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    pub directory: PathBuf,
    #[serde(default)]
    pub write_annotated_bam: bool,
}

impl Config {
    pub fn check(&self) -> Result<()> {
        self.quantification.check(self)?;
        Ok(())
    }
}

impl Input {
    pub fn read_gtf(&self) -> Result<HashMap<String, crate::gtf::GTFEntrys>> {
        if let Source::GTF(gtf_config) = &self.source {
            let accepted_features = vec![&gtf_config.feature, &gtf_config.aggr_feature];
            crate::gtf::parse_ensembl_gtf(
                &gtf_config.filename,
                accepted_features
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<HashSet<_>>(),
            )
        } else {
            bail!("Input source is not GTF, cannot read GTF entries");
        }
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
