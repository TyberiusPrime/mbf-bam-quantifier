use crate::engine::AnnotatedRead;

use super::Quant;
use anyhow::{Result, bail};
use itertools::Itertools;


use super::umi::UMIGrouping;

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct SingleCell {
    umi_grouping: UMIGrouping,
}

impl Quant for SingleCell {

    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.cell_barcodes.is_none() {
            bail!("SingleCell quantification requires cell barcodes to be defined in the configuration.");
        }
        if matches!(config.umi, crate::extractors::UMIExtraction::NoUMI(_)) {
            bail!("SingleCell quantification requires UMI extraction to be defined in the configuration.");
        }
        Ok(())
    }
    fn weight_read_group(
        &self,
        annotated_reads: &mut [(AnnotatedRead, usize)],
    ) -> Result<()> {
        //they are already sorted by (corrected) barcode.
        //All I need to do is group them and pass them off to weight_read_group
        let mut last = None;
        let mut change_indices = Vec::new();
        for (ii, (read, _org_index)) in annotated_reads.iter().enumerate() {
            match read {
                AnnotatedRead::Counted(info) => {
                    let this_barcode = Some(info.barcode.as_ref().unwrap());
                    if last.is_none() {
                        last = this_barcode;
                        change_indices.push(0);
                    } else if last != this_barcode {
                        change_indices.push(ii);
                        last = this_barcode;
                    }
                },
                _ => {unreachable!();}
            }
        }
        Ok(for (start, stop) in change_indices.iter().tuple_windows() {
            self.umi_grouping
                .weight_read_group(&mut annotated_reads[*start..*stop])
                .expect("weighting failed");
        })
    }
}


