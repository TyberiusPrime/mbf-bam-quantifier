use anyhow::{Context, Result};
use std::path::Path;

mod bam_ext;
mod barcodes;
mod categorical;
mod config;
mod deduplication;
mod engine;
mod extractors;
mod filters;
mod gtf;
mod io;
mod quantification;

use config::Config;

pub fn run(toml_file: &Path) -> Result<()> {
    let raw_config = ex::fs::read_to_string(toml_file)
        .with_context(|| format!("Could not read toml file: {}", toml_file.to_string_lossy()))?;
    let mut parsed = toml::from_str::<Config>(&raw_config)
        .with_context(|| format!("Could not parse toml file: {}", toml_file.to_string_lossy()))?;
    parsed.check().context("Error in configuration")?;
    parsed.init()?;

    quantification::quantify(
        &parsed.input,
        parsed.filter,
        &parsed.output,
        parsed.umi,
        parsed.cell_barcodes,
        parsed.dedup,
        parsed.strategy,
    )
    .context("Error in quantification")?;

    Ok(())
}
