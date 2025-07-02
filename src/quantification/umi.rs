use super::{Direction, Quant, ReadWeight};
use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Deserializer};
use std::collections::HashSet;

pub fn regex_from_string<'de, D>(deserializer: D) -> core::result::Result<regex::Regex, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    let re = regex::Regex::new(&s)
        .map_err(|e| serde::de::Error::custom(format!("Invalid regex: {e}")))?;
    Ok(re)
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedUMI {
    umi_extraction: UMIExtraction,
    umi_grouping: UMIGrouping,
    #[serde(skip)]
    storage: UMIStorage,
}

#[derive(serde::Deserialize, Debug, Clone)]
enum UMIStorage {
    Uninitialized,
    Unique(UniqueStorage),
}

impl Default for UMIStorage {
    fn default() -> Self {
        UMIStorage::Uninitialized
    }
}

impl UMIStorage {
    fn clear(&mut self) {
        match self {
            UMIStorage::Unique(storage) => {
                storage.umis_seen.clear();
            }
            UMIStorage::Uninitialized => {}
        }
    }
}

#[derive(serde::Deserialize, Debug, Clone)]
struct UniqueStorage {
    //we only need it in one directon, if read direction changes we
    //get a position_advanced
    umis_seen: HashSet<String>, //todo: replace with hashset u64, which will be enough for umis
                                //up to 32 bases...
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
enum UMIGrouping {
    #[serde(alias = "unique")]
    Unique,
}

#[enum_dispatch(UMIExtraction)]
trait UMIExtractor {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<String>;
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
#[serde(tag = "type")]
#[enum_dispatch]
enum UMIExtraction {
    #[serde(alias = "regex_name")]
    RegexName(UMIExtractionRegexName),
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
struct UMIExtractionRegexName {
    #[serde(deserialize_with = "regex_from_string")]
    regex: regex::Regex,
}

impl UMIExtractor for UMIExtractionRegexName {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<String> {
        let name = std::str::from_utf8(read.qname()).ok()?;
        self.regex
            .captures(name)
            .and_then(|caps| caps.get(1).map(|m| m.as_str().to_string()))
    }
}

impl Quant for UnstrandedUMI {
    fn init(&mut self) {
        self.storage = match self.umi_grouping {
            UMIGrouping::Unique => UMIStorage::Unique(UniqueStorage {
                umis_seen: HashSet::new(),
            }),
        }
    }

    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read(
        &mut self,
        read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        _gene_nos_seen_reverse: &HashSet<u32>,
    ) -> ReadWeight {
        let dbg_tag = read.aux(b"XF")
            .map(|v| match v {
                rust_htslib::bam::record::Aux::String(v) => v.to_string(),
                _ => panic!(),
            })
            .expect("no XF tag");
        match &mut self.storage {
            UMIStorage::Uninitialized => unreachable!(),
            UMIStorage::Unique(stor) => {
                let umi = self.umi_extraction.extract(read);
                if umi.is_none() {
                    return ReadWeight::Error(format!(
                        "umi could not be extracted from read {:?}",
                        read
                    ));
                }
                let umi = umi.unwrap();
                if dbg_tag == "ENSG00000175221.14" {
                    println!("Read {} {} has UMI: {} - was in seen: {:?}",read.pos(), std::str::from_utf8(read.qname()).unwrap(), umi,
                        stor.umis_seen.contains(&umi)
                    );
                }
                let weight = if stor.umis_seen.contains(&umi) {
                    0.
                } else {
                    stor.umis_seen.insert(umi);
                    1.0
                };

                ReadWeight::Weights(
                    gene_nos_seen_match.iter().map(|&id| (id, weight)).collect(),
                    Vec::new(),
                )
            }
        }
    }

    fn position_advanced(
        &mut self,
        _delayed_reads: &mut Vec<(rust_htslib::bam::record::Record, HashSet<u32>, HashSet<u32>)>,
        _was_final: bool,
    ) -> Vec<(
        rust_htslib::bam::record::Record,
        Vec<(u32, f64)>,
        Vec<(u32, f64)>,
    )> {
        println!("position_advanced");
        self.storage.clear();
        Vec::new()
    }
}
