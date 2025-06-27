use std::path::Path;
use anyhow::{Context, Result};

mod config;
mod bam_ext;
mod quantification;
mod categorical;
mod io;
mod gtf;

use config::Config;
use quantification::Quant;

pub fn run(toml_file: &Path) -> Result<()> {
    let raw_config = ex::fs::read_to_string(toml_file)
        .with_context(|| format!("Could not read toml file: {}", toml_file.to_string_lossy()))?;
    let mut parsed = toml::from_str::<Config>(&raw_config)
        .with_context(|| format!("Could not parse toml file: {}", toml_file.to_string_lossy()))?;
    parsed.check().context("Error in configuration")?;

    parsed
        .quantification
        .quantify(&parsed.input, &parsed.output)
        .context("Error in quantification")?;

    Ok(())
}
