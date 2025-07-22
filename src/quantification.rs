use std::collections::HashSet;

use crate::{
    config::{self, Input, Output},
    engine,
};
use anyhow::{Context, Result};
use rust_htslib::bam::Read;

fn create_output(
    mode: Option<config::OutputMode>,
    output_directory: &std::path::PathBuf,
    sorted_output_keys: Option<Vec<String>>,
    header: &str,
    first_column_only: bool,
) -> Result<engine::Output> {
    let mode = mode.expect("Should have been set by config.output.check");
    Ok(match mode {
        config::OutputMode::None => engine::Output::new_no_output(),
        config::OutputMode::Region => engine::Output::new_per_region(
            output_directory.join("counts.tsv"),
            first_column_only,
            sorted_output_keys,
            header.to_string(),
        ),
        config::OutputMode::SingleCell => {
            engine::Output::new_singlecell(output_directory.clone(), sorted_output_keys)?
        }
        config::OutputMode::StartPositions => engine::Output::new_start_positions(
            output_directory.join("start_positions.tsv"),
        ),
        config::OutputMode::Coverage => engine::Output::new_coverage(
            output_directory.join("coverage.tsv"),
        ),
    })
}

pub fn quantify(
    input: &Input,
    filters: Vec<crate::filters::Filter>,
    output: &Output,
    umi_extraction: Option<crate::config::UmiConfig>,
    cell_barcode: Option<crate::barcodes::CellBarcodes>,
    strategy: crate::config::Strategy,
) -> anyhow::Result<()> {
    // Here you would implement the quantification logic
    // For now, we just return Ok to indicate success

    // for non-cellbased counts, output only 'matching' column?
    let single_column_counts_only =
        output.only_correct || matches!(strategy.direction, crate::config::MatchDirection::Ignore);

    let our_engine = match input.source {
        crate::config::Source::Gtf(ref gtf_config) => {
            let aggr_id_attribute = gtf_config
                .aggr_id_attribute
                .as_deref()
                .unwrap_or(gtf_config.id_attribute.as_str());

            let gtf_entries = input.read_gtf(gtf_config.duplicate_handling, aggr_id_attribute)?;
            if gtf_entries.is_empty() {
                return Err(anyhow::anyhow!(
                    "No GTF entries found. Perhaps set subformat to GFF/GTF?"
                ));
            }

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
                let mut keys: Vec<String> = keys.into_iter().map(ToString::to_string).collect();
                keys.sort();
                keys
            };
            let output = create_output(
                output.mode,
                &output.directory,
                Some(sorted_output_keys.clone()),
                &aggr_id_attribute,
                single_column_counts_only,
            )?;

            engine::Engine::from_gtf(
                gtf_entries,
                gtf_config.feature.as_str(),
                gtf_config.id_attribute.as_str(),
                aggr_id_attribute,
                filters,
                umi_extraction,
                cell_barcode,
                strategy,
                output,
            )?
        }
        crate::config::Source::BamReferences => {
            let bam = rust_htslib::bam::Reader::from_path(input.bam.as_str())
                .context("Failed to open BAM file")?;
            let header = bam.header();
            let references: Result<Vec<(String, i32, u64)>> =
                header
                    .target_names()
                    .iter()
                    .enumerate()
                    .map(|(tid, name)| {
                        Ok((
                            std::str::from_utf8(name)
                                .context("reference name was'nt utf8")?
                                .to_string(),
                            tid.try_into()
                                .expect("tid did not fit in i32, but BAM spec says it's an i32"),
                            header
                                .target_len(u32::try_from(tid).expect(
                                    "tid did not fit in i32, but BAM spec says it's an i32",
                                ))
                                .context("No length for tid?!")?,
                        ))
                    })
                    .collect();
            let references = references?;
            let sorted_output_keys: Vec<String> = references
                .iter()
                .map(|(name, _tid, _length)| name.clone())
                .collect();
            let output = create_output(
                output.mode,
                &output.directory,
                Some(sorted_output_keys.clone()),
                "reference",
                single_column_counts_only,
            )?;

            engine::Engine::from_references(
                references,
                filters,
                umi_extraction,
                cell_barcode,
                strategy,
                output,
            )?
        }

        crate::config::Source::BamTag(crate::config::BamTag { tag }) => {
            let output = create_output(
                output.mode,
                &output.directory,
                None,
                std::str::from_utf8(&tag).context("Bam tag name was not valid utf8")?,
                single_column_counts_only,
            )?;

            engine::Engine::from_bam_tag(tag, filters, umi_extraction, cell_barcode, output)
        }
    };

    our_engine.quantify_bam(
        &input.bam,
        None,
        &output.directory,
        output.write_annotated_bam,
        input.max_skip_length,
        input.correct_reads_for_clipping.unwrap(),
    )?;

    Ok(())
}
