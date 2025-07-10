use itertools::izip;
use std::{
    collections::{HashMap, HashSet},
    path::PathBuf,
};

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Deserializer, Serialize};
use serde_valid::Validate;

use crate::{
    barcodes::CellBarcodes,
    deduplication::{Dedup, DeduplicationStrategy},
    extractors::UMIExtraction,
};

#[derive(Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub input: Input,
    #[serde(default)]
    pub filter: Vec<crate::filters::Filter>,
    #[serde(alias = "deduplication")]
    pub dedup: DeduplicationStrategy,

    #[serde(default)]
    pub strategy: Strategy,
    #[serde(alias = "UMI")]
    pub umi: Option<UMIExtraction>,

    pub cell_barcodes: Option<CellBarcodes>,
    pub output: Output,
}

fn default_max_skip_length() -> u32 {
    //1000 // that's what umi-tools does.
    150 // should be a decent enough value for Illumina
}

fn default_correct_reads_for_clipping() -> bool {
    true // this is the default in umi-tools
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub bam: String,
    #[serde(default = "default_correct_reads_for_clipping")]
    pub correct_reads_for_clipping: bool,
    pub source: Source,
    #[serde(default = "default_max_skip_length")]
    pub max_skip_length: u32,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
#[serde(tag = "mode")]
pub enum Source {
    #[serde(alias = "gtf")]
    Gtf(GTFConfig),
    #[serde(alias = "bam_references")]
    BamReferences,
    #[serde(alias = "bam_tag")]
    BamTag(BamTag),
}

pub fn deser_tag<'de, D>(deserializer: D) -> core::result::Result<[u8; 2], D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    let b = s.as_bytes();
    if b.len() != 2 || b[0] != b[0].to_ascii_uppercase() || b[1] != b[1].to_ascii_uppercase() {
        Err(serde::de::Error::custom(
            "Tag must be exactle two uppercase letters",
        ))?;
    }
    Ok([b[0], b[1]])
}

#[derive(Deserialize, Debug, Clone, Serialize, Validate)]
#[serde(deny_unknown_fields)]
pub struct BamTag {
    #[serde(deserialize_with = "deser_tag")]
    pub tag: [u8; 2],
}

#[derive(Deserialize, Debug, Clone, Serialize, Copy, Default)]
pub enum DuplicateHandling {
    #[serde(alias = "collapse")]
    #[default]
    Collapse,
    #[serde(alias = "rename")]
    Rename,
}

#[derive(Deserialize, Debug, Clone, Serialize, Default)]
#[serde(deny_unknown_fields)]
pub enum OverlapMode {
    #[serde(alias = "union")]
    #[default]
    Union,
    #[serde(alias = "intersection_strict")]
    IntersectionStrict,
    #[serde(alias = "intersection_non_empty")]
    IntersectionNonEmpty,
}


#[derive(Deserialize, Debug, Clone, Serialize, Default)]
#[serde(deny_unknown_fields)]
pub enum MultiRegionHandling {
    #[serde(alias = "drop")]
    Drop,
    #[serde(alias = "count_both")]
    #[default]
    CountBoth,
}


#[derive(Deserialize, Debug, Clone, Serialize, Default)]
#[serde(deny_unknown_fields)]
pub enum MatchDirection {
    #[serde(alias = "forward")]
    #[default]
    Forward,
    #[serde(alias = "reverse")]
    Reverse,
    #[serde(alias = "ignore")]
    Ignore,
}


#[derive(Deserialize, Debug, Clone, Serialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Strategy {
    #[serde(default)]
    pub overlap: OverlapMode,
    #[serde(default)]
    pub multi_region: MultiRegionHandling,
    #[serde(default)]
    pub direction: MatchDirection,
}

#[derive(Deserialize, Debug, Clone, Serialize, Copy, Default)]
#[serde(deny_unknown_fields)]
pub enum GTFFormat {
    #[default]
    AutoDetect,
    Gff,
    Gtf,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GTFConfig {
    pub filename: String,
    #[serde(default)]
    pub subformat: GTFFormat,
    pub feature: String,
    pub id_attribute: String,
    pub aggr_id_attribute: Option<String>,
    #[serde(default)]
    pub duplicate_handling: DuplicateHandling,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    pub directory: PathBuf,
    #[serde(default)]
    pub write_annotated_bam: bool,
    #[serde(default)]
    pub only_correct: bool,
}

impl Config {
    pub fn check(&self) -> Result<()> {
        self.dedup.mode.check(self)?;
        Ok(())
    }

    pub fn init(&mut self) -> Result<()> {
        if let Some(cb) = self.cell_barcodes.as_mut() {
            cb.init()?;
        }
        Ok(())
    }
}

impl Input {
    pub fn read_gtf(
        &self,
        collapse_or_rename_duplicates: crate::config::DuplicateHandling,
        duplication_detection_id_attribute: &str,
    ) -> Result<HashMap<String, crate::gtf::GTFEntrys>> {
        if let Source::Gtf(gtf_config) = &self.source {
            let accepted_features = &[&gtf_config.feature];
            let accepted_tags: HashSet<String> = vec![
                gtf_config.id_attribute.to_string(),
                gtf_config
                    .aggr_id_attribute
                    .as_ref()
                    .unwrap_or(&gtf_config.id_attribute)
                    .to_string(),
                duplication_detection_id_attribute.to_string(),
            ]
            .into_iter()
            .collect();

            /* let mut parsed = crate::gtf::parse_ensembl_gtf(
                &gtf_config.filename,
                accepted_features
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<HashSet<_>>(),
            )?; */
            let mut parsed = crate::gtf::parse_minimal(
                &gtf_config.filename,
                gtf_config.subformat,
                &accepted_features
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<HashSet<_>>(),
                &accepted_tags,
            )?;
            // gtfs tend to have repeated exons, when transcripts contain the same ones.
            // we filter those by default. But the way featureCounts does it
            // is to keep them, and then discard all their reads.
            // so for featureCount parity, we need to rename themj

            match collapse_or_rename_duplicates {
                DuplicateHandling::Collapse => {
                    for (_, entries) in parsed.iter_mut() {
                        let mut keep = vec![true; entries.seqname.len()];
                        let mut seen = HashSet::new();
                        for (ii, (start, stop, id)) in izip!(
                            &entries.start,
                            &entries.end,
                            entries
                                .vec_attributes
                                .get(duplication_detection_id_attribute)
                                .with_context(|| format!(
                                    "Attribute not found {duplication_detection_id_attribute}"
                                ))?
                        )
                        .enumerate()
                        {
                            let key = (*start, *stop, id);
                            if seen.contains(&key) {
                                keep[ii] = false;
                                continue;
                            }
                            seen.insert(key);
                        }
                        entries.filter(&keep);
                    }
                }
                DuplicateHandling::Rename => {
                    for (_, entries) in parsed.iter_mut() {
                        let mut counter = HashMap::new();
                        for (ii, (start, stop, id)) in izip!(
                            &entries.start,
                            &entries.end,
                            entries
                                .vec_attributes
                                .get_mut(duplication_detection_id_attribute)
                                .with_context(|| format!(
                                    "Attribute not found {duplication_detection_id_attribute}"
                                ))?
                        )
                        .enumerate()
                        {
                            let key = (*start, *stop, id.clone());
                            if counter.contains_key(&key) {
                                *id = format!("{}-{}", id, ii);
                                continue;
                            }
                            counter.entry(key).and_modify(|c| *c += 1).or_insert(1);
                        }
                        //entries.filter(&keep);
                    }
                }
            }

            Ok(parsed)
        } else {
            bail!("Input source is not GTF, cannot read GTF entries");
        }
    }
    /* reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        if ignore_unmapped && record.reference_sequence_id().is_none() {
            continue;
        }

        if let Some(name) = record.name() {
            func(name);
        }
    } */
}

pub fn u8_from_string<'de, D>(deserializer: D) -> core::result::Result<Vec<u8>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    Ok(s.as_bytes().to_vec())
}
