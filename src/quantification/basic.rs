use std::collections::HashSet;
use super::{Direction, Quant};

/// count every read that matches. That means a read can count for multiple
/// regions
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedBasic {}

impl Quant for UnstrandedBasic {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        let mut res = Vec::new();
        res.extend(gene_nos_seen_match.iter().map(|&id| (id, 1.0)));
        for id in gene_nos_seen_reverse {
            if !gene_nos_seen_match.contains(id) {
                res.push((*id, 1.0));
            }
        }
        (res, Vec::new())
    }
}
///
/// count every read that matches. That means a read can count for multiple
/// regions. Considers strandedness.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct StrandedBasic {
    direction: Direction,
}

impl Quant for StrandedBasic {
    fn reverse(&self) -> bool {
        match self.direction {
            Direction::Forward => false,
            Direction::Reverse => true,
        }
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        (
            gene_nos_seen_match
                .iter()
                .map(|&id| (id, 1.0))
                .collect::<Vec<_>>(),
            gene_nos_seen_reverse
                .iter()
                .map(|&id| (id, 1.0))
                .collect::<Vec<_>>(),
        )
    }
}
