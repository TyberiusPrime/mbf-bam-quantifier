use std::collections::HashSet;

use anyhow::{bail, Context, Result};
use enum_dispatch::enum_dispatch;
use rust_htslib::bam::Read;

use crate::{
    config::{Config, Input, Output},
    engine,
    extractors::UMIExtraction,
};
use serde::{Deserialize, Serialize};

mod basic;
mod singlecells;
mod umi;

#[enum_dispatch(Quantification)]
pub trait Quant: Send + Sync + Clone {
    fn check(&self, _config: &Config) -> anyhow::Result<()> {
        Ok(())
    }

    fn weight_read_group(
        &self,
        annotated_reads: &mut [(engine::AnnotatedRead, usize)],
    ) -> Result<()>;

    /// whether forward /reverse matches should be swapped
    fn reverse(&self) -> bool {
        false
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Quantification {
    #[serde(alias = "basic")]
    Basic(basic::Basic),

    #[serde(alias = "umi")]
    UMI(umi::UMI),

    #[serde(alias = "singlecell")]
    #[serde(alias = "sc")]
    SingleCell(singlecells::SingleCell),
}

impl Quantification {
    pub fn quantify(
        &mut self,
        input: &Input,
        filters: Vec<crate::filters::Filter>,
        output: &Output,
        umi_extraction: UMIExtraction,
        cell_barcode: Option<crate::barcodes::CellBarcodes>,
        strategy: crate::config::Strategy,
    ) -> anyhow::Result<()> {
        // Here you would implement the quantification logic
        // For now, we just return Ok to indicate success
        //

        let our_engine = match input.source {
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
                                        "No GTF entries found for feature {}. Perhaps set subformat to GFF/GTF? ",
                                        gtf_config.feature
                                    )
                                })?;
                    let keys: HashSet<_> = entries
                        .vec_attributes
                        .get(aggr_id_attribute)
                        .context("No aggr_id_attribute found in GTF entries")?
                        .iter()
                        .collect();
                    let mut keys: Vec<String> = keys.into_iter().map(|x| x.to_string()).collect();
                    keys.sort();
                    keys
                };

                let output = engine::Output::new_per_region(
                    output.directory.join("counts.tsv"),
                    output.only_correct || matches!(strategy.direction, crate::config::MatchDirection::Ignore),
                    Some(sorted_output_keys),
                    aggr_id_attribute.to_string(),
                );

                engine::Engine::from_gtf(
                    gtf_entries,
                    gtf_config.feature.as_str(),
                    gtf_config.id_attribute.as_str(),
                    aggr_id_attribute,
                    filters,
                    self.clone(),
                    umi_extraction,
                    cell_barcode,
                    strategy.clone(),
                    output,
                )?
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

                let output = engine::Output::new_per_region(
                    output.directory.join("counts.tsv"),
                    output.only_correct || matches!(strategy.direction, crate::config::MatchDirection::Ignore),
                    Some(sorted_output_keys),
                    "reference".to_string(),
                );

                engine::Engine::from_references(
                    references,
                    filters,
                    self.clone(),
                    umi_extraction,
                    cell_barcode,
                    strategy.clone(),
                    output,
                )?
            }

            crate::config::Source::BamTag(crate::config::BamTag { tag }) => {
                if self.reverse() {
                    bail!("Setting Direction(=reverse) on a BamTag is meaningless.")
                }

                let output = engine::Output::new_per_region(
                    output.directory.join("counts.tsv"),
                    output.only_correct || matches!(strategy.direction, crate::config::MatchDirection::Ignore),
                    None,
                    std::str::from_utf8(&tag)
                        .context("Bam tag name was not valid utf8")?
                        .to_string(),
                );
                engine::Engine::from_bam_tag(
                    tag,
                    filters,
                    self.clone(),
                    umi_extraction,
                    cell_barcode,
                    output,
                )
            }
        };

        our_engine.quantify_bam(
            &input.bam,
            None,
            &output.directory,
            output.write_annotated_bam,
            input.max_skip_length,
            input.correct_reads_for_clipping,
        )?;

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
