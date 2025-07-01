use std::{collections::HashSet, io::Write, path::Path};

use anyhow::{Context, Result};
use enum_dispatch::enum_dispatch;
use rust_htslib::bam::Read;

use crate::{
    config::{Config, Input, Output},
    engine::IntervalResult,
};
use serde::{Deserialize, Serialize};

#[enum_dispatch(Quantification)]
pub trait Quant: Send + Sync + Clone {
    fn check(&self, _config: &Config) -> anyhow::Result<()> {
        Ok(())
    }

    fn output_only_one_column(&self) -> bool {
        false
    }

    fn weight_read(
        &mut self,
        read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>);

    /// whether forward /reverse matches should be swapped
    fn reverse(&self) -> bool {
        false
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, serde::Serialize)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Quantification {
    #[serde(alias = "unstranded_basic")]
    UnstrandedBasic(UnstrandedBasic),
    #[serde(alias = "stranded_basic")]
    StrandedBasic(StrandedBasic),

    #[serde(alias = "unstranded_featurecounts")]
    UnstrandedFeatureCounts(UnstrandedFeatureCounts),
    #[serde(alias = "stranded_featurecounts")]
    StrandedFeatureCounts(StrandedFeatureCounts),
}

impl Quantification {
    pub fn quantify(
        &self,
        input: &Input,
        filters: Vec<crate::filters::Filter>,
        output: &Output,
    ) -> anyhow::Result<()> {
        // Here you would implement the quantification logic
        // For now, we just return Ok to indicate success
        //
        let (engine, sorted_output_keys, output_title) = match input.source {
            crate::config::Source::GTF(ref gtf_config) => {
                let gtf_entries = input.read_gtf(
                    gtf_config.duplicate_handling,
                    gtf_config.aggr_id_attribute.as_str(),
                )?;
                let sorted_output_keys = {
                    let entries = gtf_entries
                        .get(gtf_config.aggr_feature.as_str())
                        .with_context(|| {
                            format!("No GTF entries found for feature {}", gtf_config.feature)
                        })?;
                    let mut keys: Vec<_> = entries
                        .vec_attributes
                        .get(gtf_config.aggr_id_attribute.as_str())
                        .context("No aggr_id_attribute found in GTF entries")?
                        .iter()
                        .cloned()
                        .collect();
                    keys.sort();
                    keys
                };

                (
                    crate::engine::Engine::from_gtf(
                        gtf_entries,
                        gtf_config.feature.as_str(),
                        gtf_config.id_attribute.as_str(),
                        gtf_config.aggr_feature.as_str(),
                        gtf_config.aggr_id_attribute.as_str(),
                        filters,
                        self.clone(),
                    )?,
                    sorted_output_keys,
                    gtf_config.aggr_id_attribute.as_str(),
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
                    sorted_output_keys,
                    "reference",
                )
            }
        };

        let counts = engine.quantify_bam(
            &input.bam,
            None,
            &output.directory,
            output.write_annotated_bam,
        )?;

        Self::write_output(
            sorted_output_keys,
            counts,
            &output.directory.join("counts.tsv"),
            output_title,
            self.output_only_one_column(),
        )?;

        Ok(())
    }

    fn write_output(
        sorted_keys: Vec<String>,
        content: IntervalResult,
        output_filename: &Path,
        id_attribute: &str,
        first_column_only: bool,
    ) -> Result<()> {
        ex::fs::create_dir_all(&output_filename.parent().unwrap())?;

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

/// count every read that matches. That means a read can count for multiple
/// regions
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedBasic {}

impl Quant for UnstrandedBasic {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        let mut res = Vec::new();
        res.extend(gene_nos_seen_match.iter().map(|&id| (id, 1.0)));
        for id in gene_nos_seen_reverse {
            if !gene_nos_seen_match.contains(id) {
                res.push((*id, 1.0));
            }
        }
        (res, Vec::new())
    }
}

#[derive(Deserialize, Debug, Clone, Serialize)]
enum Direction {
    #[serde(alias = "forward")]
    Forward,
    #[serde(alias = "reverse")]
    Reverse,
}

/// count every read that matches. That means a read can count for multiple
/// regions. Considers strandedness.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct StrandedBasic {
    direction: Direction,
}

impl Quant for StrandedBasic {
    fn reverse(&self) -> bool {
        match self.direction {
            Direction::Forward => false,
            Direction::Reverse => true,
        }
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        (
            gene_nos_seen_match
                .iter()
                .map(|&id| (id, 1.0))
                .collect::<Vec<_>>(),
            gene_nos_seen_reverse
                .iter()
                .map(|&id| (id, 1.0))
                .collect::<Vec<_>>(),
        )
    }
}

/// count like featureCounts does.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedFeatureCounts {}

impl Quant for UnstrandedFeatureCounts {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        let combined = gene_nos_seen_match
            .union(gene_nos_seen_reverse)
            .cloned()
            .collect::<HashSet<_>>();
        if combined.len() == 1 {
            //matches exactly one output region.
            // Ok to count
            return (
                vec![(combined.iter().next().unwrap().clone(), 1.0)],
                Vec::new(),
            );
        } else {
            //matches multiple regions -> don't count
            (
                gene_nos_seen_match
                    .iter()
                    .map(|&id| (id, 0.0))
                    .collect::<Vec<_>>(),
                gene_nos_seen_reverse
                    .iter()
                    .map(|&id| (id, 0.0))
                    .collect::<Vec<_>>(),
            )//(Vec::new(), Vec::new())
        }
    }
}
//
// count like featureCounts does. Stranded
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct StrandedFeatureCounts {
    direction: Direction,
}

impl Quant for StrandedFeatureCounts {
    fn reverse(&self) -> bool {
        match self.direction {
            Direction::Forward => false,
            Direction::Reverse => true,
        }
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        if gene_nos_seen_match.len() == 1 {
            //matches exactly one output region.
            // Ok to count
            return (
                gene_nos_seen_match
                    .iter()
                    .map(|&id| (id, 1.0))
                    .collect::<Vec<_>>(),
                gene_nos_seen_reverse
                    .iter()
                    .map(|&id| (id, 1.0))
                    .collect::<Vec<_>>(),
            );
        } else {
            //matches multiple regions -> don't count
            (
                Vec::new(),
                gene_nos_seen_reverse
                    .iter()
                    .map(|&id| (id, 1.0))
                    .collect::<Vec<_>>(),
            )
        }
    }
}
