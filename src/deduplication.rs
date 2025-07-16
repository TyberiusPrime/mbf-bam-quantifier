use std::collections::HashMap;

use rust_htslib::bam;

use crate::{bam_ext::BamRecordExtensions};

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, Default)]
pub enum DeduplicationBucket {
    #[default]
    #[serde(alias = "per_position")]
    PerPosition,
    #[serde(alias = "per_reference")]
    PerReference,
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, Copy)]
pub enum DeduplicationMode {
    #[serde(alias = "none")]
    NoDedup,

    #[serde(alias = "unique")]
    Unique,
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

pub struct UmiGroupInfo {
    best_read_index: usize,
    best_mapping_quality: MappingQuality,
    count: usize,
}

pub enum DedupPerBucket {
    None,
    Umi(HashMap<Vec<u8>, UmiGroupInfo>),
    SingleCell(HashMap<Vec<u8>, HashMap<Vec<u8>, UmiGroupInfo>>),
}

pub enum AcceptReadResult {
    Duplicated,
    New,
    DuplicateButPrefered(usize),
}

impl DedupPerBucket {
    pub fn new(umi: bool, single_cell: bool) -> Self {
        match (umi, single_cell) {
            (true, true) => DedupPerBucket::SingleCell(HashMap::new()),
            (true, false) => DedupPerBucket::Umi(HashMap::new()),
            (false, true) => unimplemented!(), // TODO, I guess
            (false, false) => DedupPerBucket::None,
        }
    }

    pub fn accept_read(
        &mut self,
        read: &bam::record::Record,
        this_index: usize,
        umi: Option<&Vec<u8>>,
        barcode: Option<&Vec<u8>>,
        _mode: DeduplicationMode,
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
                    Some(UmiGroupInfo {
                        best_read_index: old_index,
                        best_mapping_quality: mq,
                        count,
                    }) => {
                        if this_mq > *mq {
                            *mq = this_mq;
                            let result = AcceptReadResult::DuplicateButPrefered(*old_index);
                            *old_index = this_index;
                            *count += 1;
                            result
                        } else {
                            AcceptReadResult::Duplicated
                        }
                    }
                    None => {
                        map.insert(
                            umi.to_vec(),
                            UmiGroupInfo {
                                best_read_index: this_index,
                                best_mapping_quality: this_mq,
                                count: 1,
                            },
                        );
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

                let this_mq = MappingQuality {
                    no_of_alignments: read.no_of_alignments().try_into().unwrap_or(255),
                    mapq: read.mapq(),
                };
                let by_barcode = map.entry(barcode.to_vec());
                let by_barcode = match by_barcode {
                    std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                    std::collections::hash_map::Entry::Vacant(e) => e.insert(HashMap::new()),
                };
                let hit = by_barcode.get_mut(umi);
                //todo: Combine this with above
                match hit {
                    Some(UmiGroupInfo {
                        best_read_index: old_index,
                        best_mapping_quality: mq,
                        count,
                    }) => {
                        if this_mq > *mq {
                            *mq = this_mq;
                            let result = AcceptReadResult::DuplicateButPrefered(*old_index);
                            *count += 1;
                            *old_index = this_index;
                            result
                        } else {
                            AcceptReadResult::Duplicated
                        }
                    }
                    None => {
                        by_barcode.insert(
                            umi.to_vec(),
                            UmiGroupInfo {
                                best_read_index: this_index,
                                best_mapping_quality: this_mq,
                                count: 1,
                            },
                        );
                        AcceptReadResult::New
                    }
                }
            }
        }
    }
}
