use crate::engine::AnnotatedRead;

use super::Dedup;
use anyhow::{bail, Result};
use std::collections::HashMap;

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub enum UMIGrouping {
    #[serde(alias = "unique")]
    Unique,
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct Umi {
    umi_grouping: UMIGrouping,
}

impl UMIGrouping {
    pub fn weight_read_group(&self, annotated_reads: &mut [(AnnotatedRead, usize)]) -> Result<()> {
        match self {
            UMIGrouping::Unique => {
                //just like umitools does it. group them by umi, then
                //decide on a representative read
                let mut by_umi: HashMap<&[u8], Vec<usize>> = HashMap::new();
                for (idx, (read, _org_idx)) in annotated_reads.iter().enumerate() {
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
                        //only one read with this UMI, keep it. No need ot change it.
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
                                _ => {
                                    unreachable!()
                                }
                            }
                        })
                        .unwrap());
                    for idx in indices {
                        if idx != best {
                            to_mark_duplicate.push(idx);
                        }
                    }
                }
                assert!(to_mark_duplicate.len() < annotated_reads.len());

                for idx in to_mark_duplicate {
                    annotated_reads[idx].0 = AnnotatedRead::Duplicate;
                }
            }
        }
        Ok(())
    }
}
impl Dedup for Umi {
    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.umi.is_none() {
            bail!("UMI deduplication quantification requires UMI extraction to be defined in the configuration.");
        }
        Ok(())
    }

    fn new_per_position(&self) -> super::DedupPerPosition {
        super::DedupPerPosition::Umi(HashMap::new())
    }
}
