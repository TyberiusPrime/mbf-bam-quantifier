use itertools::izip;
use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
};

use anyhow::{bail, Context, Result};
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

#[derive(Deserialize, Debug, Clone, Serialize, Copy)]
pub enum DuplicateHandling {
    #[serde(alias = "collapse")]
    Collapse,
    #[serde(alias = "rename")]
    Rename,
}

impl Default for DuplicateHandling {
    fn default() -> Self {
        DuplicateHandling::Collapse
    }
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GTFConfig {
    pub filename: String,
    pub feature: String,
    pub id_attribute: String,
    pub aggr_id_attribute: Option<String>,
    #[serde(default)]
    pub duplicate_handling: DuplicateHandling,
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
    pub fn read_gtf(
        &self,
        collapse_or_rename_duplicates: crate::config::DuplicateHandling,
        duplication_detection_id_attribute: &str,
    ) -> Result<HashMap<String, crate::gtf::GTFEntrys>> {
        if let Source::GTF(gtf_config) = &self.source {
            let accepted_features = vec![&gtf_config.feature];
            let accepted_tags: HashSet<String> = vec![
                gtf_config.id_attribute.to_string(),
                gtf_config
                    .aggr_id_attribute
                    .as_ref()
                    .unwrap_or(&gtf_config.id_attribute)
                    .to_string(),
                duplication_detection_id_attribute.to_string(),
            ]
            .into_iter()
            .collect();

            /* let mut parsed = crate::gtf::parse_ensembl_gtf(
                &gtf_config.filename,
                accepted_features
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<HashSet<_>>(),
            )?; */
            let mut parsed = crate::gtf::parse_minimal(
                &gtf_config.filename,
                accepted_features
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<HashSet<_>>(),
                accepted_tags,
            )?;
            // gtfs tend to have repeated exons, when transcripts contain the same ones.
            // we filter those by default. But the way featureCounts does it
            // is to keep them, and then discard all their reads.
            // so for featureCount parity, we need to rename themj

            match collapse_or_rename_duplicates {
                DuplicateHandling::Collapse => {
                    for (_, entries) in parsed.iter_mut() {
                        let mut keep = vec![true; entries.seqname.len()];
                        let mut seen = HashSet::new();
                        for (ii, (start, stop, id)) in izip!(
                            &entries.start,
                            &entries.end,
                            entries
                                .vec_attributes
                                .get(duplication_detection_id_attribute)
                                .with_context(|| format!(
                                    "Attribute not found {duplication_detection_id_attribute}"
                                ))?
                        )
                        .enumerate()
                        {
                            let key = (*start, *stop, id);
                            if seen.contains(&key) {
                                keep[ii] = false;
                                continue;
                            }
                            seen.insert(key);
                        }
                        entries.filter(&keep);
                    }
                }
                DuplicateHandling::Rename => {
                    for (_, entries) in parsed.iter_mut() {
                        let mut counter = HashMap::new();
                        for (ii, (start, stop, id)) in izip!(
                            &entries.start,
                            &entries.end,
                            entries
                                .vec_attributes
                                .get_mut(duplication_detection_id_attribute)
                                .with_context(|| format!(
                                    "Attribute not found {duplication_detection_id_attribute}"
                                ))?
                        )
                        .enumerate()
                        {
                            let key = (*start, *stop, id.clone());
                            if counter.contains_key(&key) {
                                *id = format!("{}-{}", id, ii);
                                continue;
                            }
                            counter.entry(key).and_modify(|c| *c += 1).or_insert(1);
                        }
                        //entries.filter(&keep);
                    }
                }
            }

            Ok(parsed)
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
