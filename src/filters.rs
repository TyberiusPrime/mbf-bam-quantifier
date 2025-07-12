use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use string_interner::symbol::SymbolU32;

use crate::bam_ext::BamRecordExtensions;
use enum_dispatch::enum_dispatch;

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize, PartialEq, Eq)]
pub enum KeepOrRemove {
    #[serde(alias = "keep")]
    Keep,
    #[serde(alias = "remove")]
    Remove,
}

#[enum_dispatch(Filter)]
pub trait ReadFilter: Send + Sync {
    fn check(&self, _config: &crate::config::Config) -> Result<()> {
        // Default implementation does nothing, can be overridden by specific filters
        Ok(())
    }

    fn init(&mut self, _header: &rust_htslib::bam::HeaderView) -> Result<()> {
        // Default implementation does nothing, can be overridden by specific filters
        Ok(())
    }

    fn remove_read(&self, _read: &rust_htslib::bam::record::Record) -> bool {
        false
    }
    fn remove_read_after_annotation(
        &self,
        _read: &rust_htslib::bam::record::Record,
        _barcode: Option<&Vec<u8>>,
        _umi: Option<&Vec<u8>>,
        _genes_hit_correct: &Vec<SymbolU32>,
        _genes_hit_reverse: &Vec<SymbolU32>,
        _interner: &crate::engine::OurInterner,
    ) -> bool {
        false
    }
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, serde::Serialize)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Filter {
    #[serde(alias = "multimapper")]
    MultiMapper(MultiMapper),
    #[serde(alias = "non_primary")]
    NonPrimary(NonPrimary),

    #[serde(alias = "read1")]
    Read1(Read1),
    #[serde(alias = "read2")]
    Read2(Read2),

    #[serde(alias = "spliced")]
    Spliced(Spliced),

    #[serde(alias = "reference")]
    Reference(Reference),

    #[serde(alias = "n_in_umi")]
    #[serde(alias = "NInUMI")]
    NInUmi(NInUmi),
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct MultiMapper {
    action: KeepOrRemove,
}

impl ReadFilter for MultiMapper {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        let alignment_count = read.no_of_alignments();
        let hit = alignment_count > 1;
        match self.action {
            KeepOrRemove::Keep => !hit,
            KeepOrRemove::Remove => hit,
        }
    }
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct NonPrimary {
    action: KeepOrRemove,
}

impl ReadFilter for NonPrimary {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        if read.is_secondary() {
            return self.action == KeepOrRemove::Remove;
        }
        false
    }
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct Read1 {
    action: KeepOrRemove,
}

impl ReadFilter for Read1 {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        if read.is_first_in_template() {
            return self.action == KeepOrRemove::Remove;
        }
        false
    }
}
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct Read2 {
    action: KeepOrRemove,
}

impl ReadFilter for Read2 {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        if read.is_last_in_template() {
            return self.action == KeepOrRemove::Remove;
        }
        false
    }
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct Spliced {
    action: KeepOrRemove,
}

impl ReadFilter for Spliced {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        let hit = read
            .cigar()
            .iter()
            .skip(1)
            .any(|c| matches!(c, rust_htslib::bam::record::Cigar::RefSkip(_)));
        match self.action {
            KeepOrRemove::Keep => !hit,
            KeepOrRemove::Remove => hit,
        }
    }
}

/// Filter references from processing.
/// Note that unlike the other filters, this one will even filter them from output
/// since we're leveraging this to filter the chunks we're looking at.
/// (This way, we're both saving the time and the memory to process them,
/// which can be a life safer if you've a file with a gazillion reads mapping to the same
/// coordinates, for example phiX).
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct Reference {
    pub action: KeepOrRemove,
    pub references: Vec<String>,
    #[serde(skip)]
    tids: Option<HashSet<u32>>,
}

impl ReadFilter for Reference {
    fn check(&self, _config: &crate::config::Config) -> Result<()> {
        if self.references.is_empty() {
            bail!("Reference filter requires at least one reference to filter on.");
        }
        Ok(())
    }

    fn init(&mut self, header: &rust_htslib::bam::HeaderView) -> Result<()> {
        let mut tids = HashSet::new();
        for r in &self.references {
            let tid = header
                .tid(r.as_bytes())
                .with_context(|| format!("Reference not found in header: {r}"))?;
            tids.insert(tid);
        }
        self.tids = Some(tids);
        Ok(())
    }
    fn remove_read(&self, _read: &rust_htslib::bam::record::Record) -> bool {
        false
        /* has already been filtered by the chunk filtering
                let hit = self.tids.as_ref().unwrap().contains(&(read.tid() as u32));
                match self.action {
                    KeepOrRemove::Keep => !hit,
                    KeepOrRemove::Remove => hit,
                }
        */
    }
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct NInUmi {
    pub action: KeepOrRemove,
}

impl ReadFilter for NInUmi {
    fn remove_read_after_annotation(
        &self,
        _read: &rust_htslib::bam::record::Record,
        _barcode: Option<&Vec<u8>>,
        umi: Option<&Vec<u8>>,
        _genes_hit_correct: &Vec<SymbolU32>,
        _genes_hit_reverse: &Vec<SymbolU32>,
        _interner: &crate::engine::OurInterner,
    ) -> bool {
        if let Some(umi) = umi {
            let hit = umi.iter().any(|x| *x == b'N');

            match self.action {
                KeepOrRemove::Keep => !hit,
                KeepOrRemove::Remove => hit,
            }
        } else {
            false
        }
    }
}
