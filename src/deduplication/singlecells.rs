use std::collections::HashMap;


use super::Dedup;
use anyhow::{bail};

use super::umi::UMIGrouping;

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct SingleCell {
    umi_grouping: UMIGrouping,
}

impl Dedup for SingleCell {
    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.cell_barcodes.is_none() {
            bail!("SingleCell quantification requires cell barcodes to be defined in the configuration.");
        }
        if config.umi.is_none() {
            bail!("SingleCell quantification requires UMI extraction to be defined in the configuration.");
        }
        Ok(())
    }

}
