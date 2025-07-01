use std::{
    collections::{HashMap, HashSet},
    io::Write,
    path::Path,
};

use anyhow::{bail, Context, Result};
use enum_dispatch::enum_dispatch;

use crate::{
    config::{Config, Input, Output},
    engine,
    engine::IntervalResult,
};

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
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, serde::Serialize)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Quantification {
    #[serde(alias = "unstranded_basic")]
    UnstrandedBasic(UnstrandedBasic),
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
        let gtf_entries = input.read_gtf()?;
        let gtf_config = input.gtf.as_ref().context({
            anyhow::anyhow!("GTF file is required for quantification, but none was provided.")
        })?;
        let mut sorted_output_keys = {
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

        let engine = crate::engine::Engine::from_gtf(
            gtf_entries,
            gtf_config.feature.as_str(),
            gtf_config.id_attribute.as_str(),
            gtf_config.aggr_feature.as_str(),
            gtf_config.aggr_id_attribute.as_str(),
            filters,
            self.clone(),
        )?;

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
            gtf_config.aggr_id_attribute.as_str(),
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

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct UnstrandedBasic {}

impl Quant for UnstrandedBasic {
    fn check(&self, config: &Config) -> anyhow::Result<()> {
        if config.input.gtf.is_none() {
            bail!("UnstrandedBasic quantification requires a GTF file in the input configuration.");
        }
        Ok(())
    }

    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read(
        &mut self,
        read: &rust_htslib::bam::record::Record,
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
