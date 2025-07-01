use anyhow::{Context, Result};
use std::{path::Path};

mod bam_ext;
mod categorical;
mod config;
mod engine;
mod filters;
mod gtf;
mod io;
mod quantification;

use config::Config;

pub fn run(toml_file: &Path) -> Result<()> {
    let raw_config = ex::fs::read_to_string(toml_file)
        .with_context(|| format!("Could not read toml file: {}", toml_file.to_string_lossy()))?;
    let parsed = toml::from_str::<Config>(&raw_config)
        .with_context(|| format!("Could not parse toml file: {}", toml_file.to_string_lossy()))?;
    parsed.check().context("Error in configuration")?;

    parsed
        .quantification
        .quantify(&parsed.input, parsed.filter, &parsed.output)
        .context("Error in quantification")?;

    Ok(())
}
