use crate::engine::AnnotatedRead;

use super::{Direction, Quant};
use anyhow::Result;
use serde::{Deserialize, Deserializer};
use std::collections::{HashMap, HashSet};



#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedUMI {
    umi_grouping: UMIGrouping,
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
enum UMIGrouping {
    #[serde(alias = "unique")]
    Unique,
}

impl Quant for UnstrandedUMI {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read_group(
        &self,
        annotated_reads: &mut [(crate::engine::AnnotatedRead, usize)],
    ) -> Result<()> {
        //I have a read group. They all start at the same (corrected) position.
        match self.umi_grouping {
            UMIGrouping::Unique => {
                //just like umitools does it. group them by umi, then
                //decide on a representative read
                let mut by_umi: HashMap<&str, Vec<usize>> = HashMap::new();
                for (idx, (read, org_idx)) in annotated_reads.iter().enumerate() {
                    if let AnnotatedRead::Counted(info) = read {
                        //extract the UMI
                        let umi = info
                            .umi
                            .as_ref()
                            .expect("UMIExtractor should have extracted UMI");
                        by_umi.entry(umi).or_default().push(idx);
                    }
                }

                let mut to_mark_duplicate = Vec::new();
                for (_, indices) in by_umi {
                    if indices.len() == 1 {
                        //only one read with this UMI, keep it. No need ot achange it.
                        continue;
                    }
                    // multiple reads with the same UMI
                    // we follow almost the umi-tools schema
                    // but no randomness!
                    let best = *(indices
                        .iter()
                        .max_by_key(|idx| {
                            let (anno_read, org_index) = &annotated_reads[**idx];
                            match anno_read {
                                AnnotatedRead::Counted(info) => (info.mapping_priority, org_index),
                                _ => {unreachable!()}
                            }
                        })
                        .unwrap());
                    for idx in indices {
                        if idx != best {
                            to_mark_duplicate.push(idx);
                        }
                    }
                }

                for idx in to_mark_duplicate {
                    annotated_reads[idx].0 = AnnotatedRead::Duplicate;
                }


            }
        }
        Ok(())
    }
}
