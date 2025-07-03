use enum_dispatch::enum_dispatch;
use crate::bam_ext::BamRecordExtensions;

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize, PartialEq,Eq)]
enum KeepOrRemove {
    #[serde(alias = "keep")]
    Keep,
    #[serde(alias = "remove")]
    Remove,
}

#[enum_dispatch(Filter)]
pub trait ReadFilter: Send + Sync {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool;
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
