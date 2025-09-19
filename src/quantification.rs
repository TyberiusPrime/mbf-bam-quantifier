use std::collections::{HashMap, HashSet};

use crate::{
    config::{self, Input, Output},
    engine,
};
use anyhow::{bail, Context, Result};
use ex::Wrapper;
use itertools::izip;
use rust_htslib::bam::Read;
use std::io::Write;

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
        config::OutputMode::StartPositions => {
            engine::Output::new_start_positions(output_directory.join("start_positions.tsv"))
        }
        config::OutputMode::Coverage => {
            engine::Output::new_coverage(output_directory.join("coverage.tsv"))
        }
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
            let engine_output = create_output(
                output.mode,
                &output.directory,
                Some(sorted_output_keys.clone()),
                &aggr_id_attribute,
                single_column_counts_only,
            )?;

            if output.output_effective_lengths {
                let lens  = estimate_effective_lengths(
                            gtf_entries
                                .get(gtf_config.feature.as_str())
                                .with_context(|| {
                                    format!(
                                        "No GTF entries found for feature {}. Perhaps set subformat to GFF/GTF? ",
                                        gtf_config.feature
                                    )
                                })?,
                aggr_id_attribute,
                )?;
                write_lengths(
                    &lens,
                    aggr_id_attribute,
                    &sorted_output_keys,
                    &output.directory.join("effective_lengths.tsv"),
                )?;
            }

            let e = engine::Engine::from_gtf(
                gtf_entries,
                gtf_config.feature.as_str(),
                gtf_config.id_attribute.as_str(),
                aggr_id_attribute,
                filters,
                umi_extraction,
                cell_barcode,
                strategy,
                engine_output,
            )?;
            e
        }
        crate::config::Source::BamReferences => {
            if output.output_effective_lengths {
                bail!("TODO: output effective lengths when using BAM references");
            }
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

            if output.output_effective_lengths {
                bail!(
                    "TODO: can not output effective lengths when using BAM references. Turn off output.output_effective_lengths"
                );
            }
            let engine_output = create_output(
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
                engine_output,
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

fn write_lengths(
    lens: &HashMap<String, usize>,
    aggregation_id_attribute: &str,
    sorted_output_keys: &[String],
    output_filename: &std::path::PathBuf,
) -> Result<()> {
    let mut handle = std::io::BufWriter::new(ex::fs::File::create(output_filename)?.into_inner());
    handle
        .write_all(format!("{aggregation_id_attribute}\tlength\n").as_bytes())
        .context("Failed to write to file")?;

    for key in sorted_output_keys {
        let len = lens.get(key).context("Missing entry in lengths, but in sorted output keys???")?;
        handle
            .write_all(format!("{key}\t{len}\n").as_bytes())
            .context("Failed to write to file")?;
    }

    Ok(())

}

fn estimate_effective_lengths(
    gtf_entries: &crate::gtf::GTFEntrys,
    aggregation_id_attribute: &str,
) -> Result<HashMap<String, usize>>{
    let mut intervals = HashMap::new();
    for (_seq_name_cat_id, gene_id, start, end, _strand) in izip!(
        gtf_entries.seqname.values.iter(),
        gtf_entries
            .vec_attributes
            .get(aggregation_id_attribute)
            .context("Missing aggr_id attribute")?
            .iter(),
        gtf_entries.start.iter(),
        gtf_entries.end.iter(),
        gtf_entries.strand.iter(),
    ) {
        let key = gene_id.clone();
        let entry = intervals.entry(key).or_insert_with(Vec::new);
        entry.push(*start as u32..*end as u32);
    }
    let result = intervals
        .into_iter()
        .map(|(gene_id, intervals)| {
            let ivset = nested_intervals::IntervalSet::new(&intervals).unwrap();
            (gene_id, ivset.covered_units() as usize)
        })
        .collect();
    Ok(result)
}
