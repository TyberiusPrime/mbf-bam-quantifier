use anyhow::Result;
use enum_dispatch::enum_dispatch;

use crate::{config::Config, engine};
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
