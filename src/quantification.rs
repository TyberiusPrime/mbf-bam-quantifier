use std::{collections::HashSet, io::Write, path::Path};

use anyhow::{bail, Context, Result};
use enum_dispatch::enum_dispatch;
use rust_htslib::bam::Read;

use crate::{
    config::{Config, Input, Output},
    engine::IntervalResult,
};
use serde::{Deserialize, Serialize};

mod basic;
mod featurecounts;
/*mod umi; */

#[enum_dispatch(Quantification)]
pub trait Quant: Send + Sync + Clone {
    fn check(&self, _config: &Config) -> anyhow::Result<()> {
        Ok(())
    }
    fn init(&mut self) {}

    fn output_only_one_column(&self) -> bool {
        false
    }

    fn weight_read_group(&self, annotated_reads: &mut [crate::engine::AnnotatedRead])
        -> Result<()>;

    /// whether forward /reverse matches should be swapped
    fn reverse(&self) -> bool {
        false
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Quantification {
    #[serde(alias = "unstranded_basic")]
    UnstrandedBasic(basic::UnstrandedBasic),

    #[serde(alias = "stranded_basic")]
    StrandedBasic(basic::StrandedBasic),

    #[serde(alias = "unstranded_featurecounts")]
    UnstrandedFeatureCounts(featurecounts::UnstrandedFeatureCounts),
    #[serde(alias = "stranded_featurecounts")]
    StrandedFeatureCounts(featurecounts::StrandedFeatureCounts),
    /*

    #[serde(alias = "unstranded_umi")]
    UnstrandedUMI(umi::UnstrandedUMI), */
}

impl Quantification {
    pub fn quantify(
        &mut self,
        input: &Input,
        filters: Vec<crate::filters::Filter>,
        output: &Output,
    ) -> anyhow::Result<()> {
        // Here you would implement the quantification logic
        // For now, we just return Ok to indicate success
        //
        self.init();
        let (engine, sorted_output_keys, output_title): (_, Option<Vec<String>>, String) =
            match input.source {
                crate::config::Source::GTF(ref gtf_config) => {
                    let aggr_id_attribute = gtf_config
                        .aggr_id_attribute
                        .as_ref()
                        .map(|x| x.as_str())
                        .unwrap_or(gtf_config.id_attribute.as_str());

                    let gtf_entries =
                        input.read_gtf(gtf_config.duplicate_handling, aggr_id_attribute)?;

                    let sorted_output_keys = {
                        let entries =
                            gtf_entries
                                .get(gtf_config.feature.as_str())
                                .with_context(|| {
                                    format!(
                                        "No GTF entries found for feature {}",
                                        gtf_config.feature
                                    )
                                })?;
                        let keys: HashSet<_> = entries
                            .vec_attributes
                            .get(aggr_id_attribute)
                            .context("No aggr_id_attribute found in GTF entries")?
                            .iter()
                            .collect();
                        let mut keys: Vec<String> =
                            keys.into_iter().map(|x| x.to_string()).collect();
                        keys.sort();
                        keys
                    };

                    (
                        crate::engine::Engine::from_gtf(
                            gtf_entries,
                            gtf_config.feature.as_str(),
                            gtf_config.id_attribute.as_str(),
                            aggr_id_attribute,
                            filters,
                            self.clone(),
                        )?,
                        Some(sorted_output_keys),
                        aggr_id_attribute.to_string(),
                    )
                }
                crate::config::Source::BamReferences => {
                    let bam = rust_htslib::bam::Reader::from_path(input.bam.as_str())
                        .context("Failed to open BAM file")?;
                    let header = bam.header();
                    let references: Result<Vec<(String, u64)>> = header
                        .target_names()
                        .iter()
                        .enumerate()
                        .map(|(tid, name)| {
                            Ok((
                                std::str::from_utf8(name)
                                    .context("reference name was'nt utf8")?
                                    .to_string(),
                                header
                                    .target_len(tid as u32)
                                    .context("No length for tid?!")?,
                            ))
                        })
                        .collect();
                    let references = references?;
                    let sorted_output_keys: Vec<String> =
                        references.iter().map(|(name, _)| name.clone()).collect();
                    (
                        crate::engine::Engine::from_references(references, filters, self.clone())?,
                        Some(sorted_output_keys),
                        "reference".to_string(),
                    )
                }

                crate::config::Source::BamTag(crate::config::BamTag { tag }) => {
                    if self.reverse() {
                        bail!("Setting Direction(=reverse) on a BamTag is meaningless.")
                    }
                    (
                        crate::engine::Engine::from_bam_tag(tag, filters, self.clone()),
                        None,
                        std::str::from_utf8(&tag)
                            .context("Bam tag was not valid utf8")?
                            .to_string(),
                    )
                }
            };

        let counts = engine.quantify_bam(
            &input.bam,
            None,
            &output.directory,
            output.write_annotated_bam,
            input.max_skip_length,
        )?;

        Self::write_output(
            sorted_output_keys,
            counts,
            &output.directory.join("counts.tsv"),
            &output_title,
            self.output_only_one_column(),
        )?;

        Ok(())
    }

    fn write_output(
        sorted_keys: Option<Vec<String>>,
        content: IntervalResult,
        output_filename: &Path,
        id_attribute: &str,
        first_column_only: bool,
    ) -> Result<()> {
        ex::fs::create_dir_all(&output_filename.parent().unwrap())?;

        let sorted_keys = sorted_keys.unwrap_or_else(|| {
            let mut keys: Vec<_> = content.counts.keys().cloned().collect();
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
                let count = content.counts.get(&key).unwrap_or(&(0.0, 0.0)).0;
                out_buffer
                    .write_all(format!("{}\t{}\n", key, count).as_bytes())
                    .context("Failed to write counts to output file")?;
            }
        } else {
            out_buffer
                .write_all(format!("{}\tcount_correct\tcount_reverse\n", id_attribute).as_bytes())
                .context("Failed to write header to output file")?;

            for key in sorted_keys {
                let (count_correct, count_reverse) =
                    content.counts.get(&key).unwrap_or(&(0.0, 0.0));
                out_buffer
                    .write_all(
                        format!("{}\t{}\t{}\n", key, count_correct, count_reverse).as_bytes(),
                    )
                    .context("Failed to write counts to output file")?;
            }
        }

        Ok(())
    }
}

#[derive(Deserialize, Debug, Clone, Serialize)]
enum Direction {
    #[serde(alias = "forward")]
    Forward,
    #[serde(alias = "reverse")]
    Reverse,
}
