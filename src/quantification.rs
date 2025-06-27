mod contig;

use enum_dispatch::enum_dispatch;
use crate::config::{Input, Output, Config};

#[enum_dispatch(Quantification)]
pub trait Quant {
    fn quantify(&mut self, input: &Input, output: &Output) -> anyhow::Result<()>;
    fn check(&self, config: &Config) -> anyhow::Result<()> {Ok(())}
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, serde::Serialize)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Quantification {
    Contig(contig::Quantification),
}

 
impl Quantification {
}
