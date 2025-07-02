/// Quantifications that work like featureCounts
/// from the subread package does.
use super::{Direction, Quant};
use std::collections::HashSet;

/// count like featureCounts does.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedFeatureCounts {}

impl Quant for UnstrandedFeatureCounts {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read(
        &mut self,
        _read: &rust_htslib::bam::record::Record,
        gene_nos_seen_match: &HashSet<u32>,
        gene_nos_seen_reverse: &HashSet<u32>,
    ) -> (Vec<(u32, f64)>, Vec<(u32, f64)>) {
        let combined = gene_nos_seen_match
            .union(gene_nos_seen_reverse)
            .cloned()
            .collect::<HashSet<_>>();
        if combined.len() == 1 {
            //matches exactly one output region.
            // Ok to count
            return (
                vec![(combined.iter().next().unwrap().clone(), 1.0)],
                Vec::new(),
            );
        } else {
            //matches multiple regions -> don't count
            (
                gene_nos_seen_match
                    .iter()
                    .map(|&id| (id, 0.0))
                    .collect::<Vec<_>>(),
                gene_nos_seen_reverse
                    .iter()
                    .map(|&id| (id, 0.0))
                    .collect::<Vec<_>>(),
            )//(Vec::new(), Vec::new())
        }
    }
}

//
// count like featureCounts does. Stranded
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct StrandedFeatureCounts {
    direction: Direction,
}

impl Quant for StrandedFeatureCounts {
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
        if gene_nos_seen_match.len() == 1 {
            //matches exactly one output region.
            // Ok to count
            return (
                gene_nos_seen_match
                    .iter()
                    .map(|&id| (id, 1.0))
                    .collect::<Vec<_>>(),
                gene_nos_seen_reverse
                    .iter()
                    .map(|&id| (id, 1.0))
                    .collect::<Vec<_>>(),
            );
        } else {
            //matches multiple regions -> don't count
            (
                Vec::new(),
                gene_nos_seen_reverse
                    .iter()
                    .map(|&id| (id, 1.0))
                    .collect::<Vec<_>>(),
            )
        }
    }
}
