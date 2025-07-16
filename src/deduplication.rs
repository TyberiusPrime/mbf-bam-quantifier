use std::collections::{HashMap, HashSet};

use rust_htslib::bam;

use crate::{bam_ext::BamRecordExtensions, engine};

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

    #[serde(alias = "percentile")]
    Percentile,
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

    pub fn accept_read_inner(
        map: &mut HashMap<Vec<u8>, UmiGroupInfo>,
        this_index: usize,
        read: &bam::record::Record,
        umi: &Vec<u8>,
    ) -> AcceptReadResult {
        let hit = map.get_mut(umi);
        let this_mq = MappingQuality {
            no_of_alignments: read.no_of_alignments().try_into().unwrap_or(255),
            mapq: read.mapq(),
        };
        match hit {
            Some(UmiGroupInfo {
                best_read_index: old_index,
                best_mapping_quality: mq,
                count,
            }) => {
                *count += 1;
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

    pub fn accept_read(
        &mut self,
        read: &bam::record::Record,
        this_index: usize,
        umi: Option<&Vec<u8>>,
        barcode: Option<&Vec<u8>>,
    ) -> AcceptReadResult {
        match self {
            DedupPerBucket::None => AcceptReadResult::New,
            DedupPerBucket::Umi(map) => Self::accept_read_inner(
                map,
                this_index,
                read,
                umi.expect("UMI should be extracted before deduplication"),
            ),
            DedupPerBucket::SingleCell(map) => {
                let barcode = barcode
                    .expect("Barcode should be extracted before deduplication")
                    .as_slice();

                let by_barcode = map.entry(barcode.to_vec());
                let by_barcode = match by_barcode {
                    std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                    std::collections::hash_map::Entry::Vacant(e) => e.insert(HashMap::new()),
                };
                Self::accept_read_inner(
                    by_barcode,
                    this_index,
                    read,
                    umi.expect("UMI should be extracted before deduplication"),
                )
            }
        }
    }

    pub fn finish_bucket(
        &mut self,
        reads: &mut Vec<(engine::AnnotatedRead, usize)>,
        mode: DeduplicationMode,
    ) {
        match self {
            DedupPerBucket::None => {}
            DedupPerBucket::Umi(map) => {
                Self::finish_umi_subbucket(mode, &map, reads);
            }
            DedupPerBucket::SingleCell(by_barcode) => {
                for (_barcode, map) in by_barcode.iter_mut() {
                    Self::finish_umi_subbucket(mode, map, reads);
                }
            }
        }
    }

    fn finish_umi_subbucket(
        mode: DeduplicationMode,
        map: &HashMap<Vec<u8>, UmiGroupInfo>,
        reads: &mut Vec<(engine::AnnotatedRead, usize)>,
    ) {
        match mode {
            DeduplicationMode::NoDedup => {}
            DeduplicationMode::Unique => {}
            DeduplicationMode::Percentile => {
                //UMI. UMIs with counts < 1% of the median counts for UMIs at the same position are ignored.
                const _: () = assert!(
                    std::mem::size_of::<usize>() == 8,
                    "usize isn't u64, this needs adjusting"
                );
                let mut counts: Vec<u64> = map
                    .values()
                    .map(
                        |v| v.count as u64, //if usize > u64
                        //the const assert above should have exploded
                    )
                    .collect();
                if counts.len() == 1 {
                    return
                }
                let medians =
                    medians::medianu64(&mut counts[..]).unwrap_or(medians::Medians::Odd(&0));
                let median = match medians {
                    medians::Medians::Odd(m) => *m,
                    medians::Medians::Even((m1, m2)) => (m1 + m2) / 2,
                };

                let threshold = median / 100; //rounds towards 0, so should be fine.
                if counts.iter().any(|x| *x > 1) {
                    dbg!(&counts);
                    dbg!(median);
                    dbg!(&threshold);
                }
                let filter: HashSet<_> = map
                    .iter()
                    .filter(|(_k, v)| (v.count as u64) <= threshold)
                    .map(|(k, _v)| k)
                    .collect();
                for (read, _org_index) in reads.iter_mut() {
                    match read {
                        engine::AnnotatedRead::Counted(info) => {
                            if let Some(umi) = info.umi.as_ref() {
                                if filter.contains(umi) {
                                    *read = engine::AnnotatedRead::new_approximate_duplicate(
                                        info.corrected_position,
                                    );
                                }
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
    }
}
