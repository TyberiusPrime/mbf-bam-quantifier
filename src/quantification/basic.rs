use crate::engine::AnnotatedRead;

use super::{Direction, Quant};
use anyhow::Result;

/// count every read that matches. That means a read can count for multiple
/// regions
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedBasic {}

impl Quant for UnstrandedBasic {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read_group(
        &self,
        annotated_reads: &mut [(crate::engine::AnnotatedRead, usize)],
    ) -> Result<()> {
        for (read, org_index) in annotated_reads.iter_mut() {
            match read {
                AnnotatedRead::Counted(read) => {
                    for gene_id in read.genes_hit_reverse.keys() {
                        if !read.genes_hit_correct.contains_key(gene_id) {
                            read.genes_hit_correct.insert(gene_id.to_string(), 1.0);
                        }
                    }
                    read.genes_hit_reverse.clear()
                }
                _ => {
                    unreachable!();
                }
            }
        }
        Ok(())
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
    fn weight_read_group(
        &self,
        _annotated_reads: &mut [(crate::engine::AnnotatedRead, usize)],
    ) -> Result<()> {
        //since we already have 1.0 in the forward direction
        //for simple matching reads, this is a noop.
        //
        Ok(())
    }
}
