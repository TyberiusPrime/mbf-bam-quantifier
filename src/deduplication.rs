use std::collections::HashMap;

use anyhow::Result;
use enum_dispatch::enum_dispatch;
use rust_htslib::bam;

use crate::{bam_ext::BamRecordExtensions, config::Config, engine};
use serde::{Deserialize, Serialize};

mod basic;
mod singlecells;
mod umi;

#[enum_dispatch(DeduplicationStrategy)]
pub trait Dedup: Send + Sync + Clone {
    fn check(&self, _config: &Config) -> anyhow::Result<()> {
        Ok(())
    }

    fn dedup(&self, annotated_reads: &mut [(engine::AnnotatedRead, usize)]) -> Result<()>;

    fn new_per_position(&self) -> DedupPerPosition {
        DedupPerPosition::None
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum DeduplicationStrategy {
    #[serde(alias = "none")]
    NoDedup(basic::Basic), //todo: rename basic:basic

    #[serde(alias = "umi")]
    UMI(umi::UMI),

    #[serde(alias = "singlecell")]
    #[serde(alias = "sc")]
    SingleCell(singlecells::SingleCell),
}

#[derive(Deserialize, Debug, Clone, Serialize)]
enum Direction {
    #[serde(alias = "forward")]
    Forward,
    #[serde(alias = "reverse")]
    Reverse,
}

#[derive(PartialOrd, PartialEq, Eq)]
struct MappingQuality {
    no_of_alignments: u8,
    mapq: u8,
    //consider deciding on cigar length as well?
}

impl Ord for MappingQuality {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.no_of_alignments
            .cmp(&other.no_of_alignments)
            .then(self.mapq.cmp(&other.mapq))
    }
}

pub enum DedupPerPosition {
    None,
    UMI(HashMap<Vec<u8>, (usize, MappingQuality)>),
    SingleCell(HashMap<(Vec<u8>, Vec<u8>), (usize, MappingQuality)>),
}

pub enum AcceptReadResult {
    Duplicated,
    New,
    DuplicateButPrefered(usize),
}

impl DedupPerPosition {
    pub fn accept_read(
        &mut self,
        read: &bam::record::Record,
        this_index: usize,
        umi: Option<&Vec<u8>>,
        barcode: Option<&Vec<u8>>,
    ) -> AcceptReadResult {
        match self {
            DedupPerPosition::None => return AcceptReadResult::New,
            DedupPerPosition::UMI(map) => {
                let umi = umi
                    .expect("UMI should be extracted before deduplication")
                    .as_slice();
                let this_mq = MappingQuality {
                    no_of_alignments: read.no_of_alignments().try_into().unwrap_or(255),
                    mapq: read.mapq(),
                };
                let hit = map.get_mut(umi);
                match hit {
                    Some((old_index, mq)) => {
                        if this_mq > *mq {
                            *mq = this_mq;
                            let result = AcceptReadResult::DuplicateButPrefered(*old_index);
                            *old_index = this_index;
                            return result;
                        } else {
                            AcceptReadResult::Duplicated
                        }
                    }
                    None => {
                        map.insert(umi.to_vec(), (this_index, this_mq));
                        return AcceptReadResult::New;
                    }
                }
            }
            DedupPerPosition::SingleCell(map) => {
                let umi = umi
                    .expect("UMI should be extracted before deduplication")
                    .as_slice();
                let barcode = barcode
                    .expect("Barcode should be extracted before deduplication")
                    .as_slice();

                let this_mq = MappingQuality {
                    no_of_alignments: read.no_of_alignments().try_into().unwrap_or(255),
                    mapq: read.mapq(),
                };
                let hit = map.get_mut(&(umi.to_vec(), barcode.to_vec()));
                match hit {
                    Some((old_index, mq)) => {
                        if this_mq > *mq {
                            *mq = this_mq;
                            let result = AcceptReadResult::DuplicateButPrefered(*old_index);
                            *old_index = this_index;
                            return result;
                        } else {
                            AcceptReadResult::Duplicated
                        }
                    }
                    None => 
                        map.insert((umi.to_vec(), barcode.to_vec()), (this_index, this_mq));
                        return AcceptReadResult::New;
                    }
                }
            }
        }
    }
}
