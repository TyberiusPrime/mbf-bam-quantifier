use std::collections::HashMap;

use anyhow::bail;
use rust_htslib::bam;

use crate::{bam_ext::BamRecordExtensions, config::Config};

#[derive(serde::Deserialize, Debug, Clone)]
pub struct DeduplicationStrategy {
    //this is 'flattend' into the [dedup]
    #[serde(flatten)]
    pub mode: DeduplicationMode,
    #[serde(default)]
    pub bucket: DeduplicationBucket,
}

impl DeduplicationStrategy {
    pub fn new_bucket(&self) -> DedupPerBucket {
        match &self.mode {
            DeduplicationMode::NoDedup => DedupPerBucket::None,
            DeduplicationMode::Umi => DedupPerBucket::Umi(HashMap::new()),
            DeduplicationMode::SingleCell => DedupPerBucket::SingleCell(HashMap::new()),
        }
    }

    pub fn check(&self, config: &Config) -> anyhow::Result<()> {
        match self.mode {
            DeduplicationMode::NoDedup => {}
            DeduplicationMode::Umi => {
                if config.umi.is_none() {
                    bail!("UMI deduplication quantification requires UMI extraction to be defined in the configuration.");
                }
            }
            DeduplicationMode::SingleCell => {
                if config.cell_barcodes.is_none() {
                    bail!("SingleCell quantification requires cell barcodes to be defined in the configuration.");
                }
                if config.umi.is_none() {
                    bail!("SingleCell quantification requires UMI extraction to be defined in the configuration.");
                }
            }
        }
        Ok(())
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, Default)]
pub enum DeduplicationBucket {
    #[default]
    PerPosition,
    PerReference,
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display)]
//have it read the mode field
#[serde(tag = "mode")]
pub enum DeduplicationMode {
    #[serde(alias = "none")]
    NoDedup,

    #[serde(alias = "umi")]
    Umi,

    #[serde(alias = "singlecell")]
    #[serde(alias = "sc")]
    SingleCell,
}

#[derive(PartialEq, Eq, Debug)]
pub struct MappingQuality {
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

impl PartialOrd for MappingQuality {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub enum DedupPerBucket {
    None,
    Umi(HashMap<Vec<u8>, (usize, MappingQuality)>),
    SingleCell(HashMap<(Vec<u8>, Vec<u8>), (usize, MappingQuality)>),
}

pub enum AcceptReadResult {
    Duplicated,
    New,
    DuplicateButPrefered(usize),
}

impl DedupPerBucket {
    pub fn accept_read(
        &mut self,
        read: &bam::record::Record,
        this_index: usize,
        umi: Option<&Vec<u8>>,
        barcode: Option<&Vec<u8>>,
    ) -> AcceptReadResult {
        match self {
            DedupPerBucket::None => AcceptReadResult::New,
            DedupPerBucket::Umi(map) => {
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
                            result
                        } else {
                            AcceptReadResult::Duplicated
                        }
                    }
                    None => {
                        map.insert(umi.to_vec(), (this_index, this_mq));
                        AcceptReadResult::New
                    }
                }
            }
            DedupPerBucket::SingleCell(map) => {
                let umi = umi
                    .expect("UMI should be extracted before deduplication")
                    .as_slice();
                let barcode = barcode
                    .expect("Barcode should be extracted before deduplication")
                    .as_slice();

                /* println!("read {:?} pos {}, umi {}, barcode {}",
                                   std::str::from_utf8(read.qname()),
                                   read.pos(),
                                   std::str::from_utf8(umi).unwrap_or("invalid UMI"),
                                   std::str::from_utf8(barcode).unwrap_or("invalid barcode"));
                */
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
                            result
                        } else {
                            AcceptReadResult::Duplicated
                        }
                    }
                    None => {
                        map.insert((umi.to_vec(), barcode.to_vec()), (this_index, this_mq));
                        AcceptReadResult::New
                    }
                }
            }
        }
    }
}
