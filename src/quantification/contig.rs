use super::Quant;
use crate::config::{Input, Output};
use anyhow::bail;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Quantification {}

impl Quant for Quantification {
    fn quantify(&mut self, input: &Input, output: &Output) -> anyhow::Result<()> {
        // Implement the quantification logic here
        todo!();
        Ok(())
    }

    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.input.gtf.is_some() {
            bail!("GTF file is not required for contig quantification, do net set it.");
        }

        Ok(())
    }
}
