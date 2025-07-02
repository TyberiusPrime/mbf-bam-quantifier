/// Quantifications that work like featureCounts
/// from the subread package does.
use super::{Direction, Quant};
use anyhow::Result;

/// count like featureCounts does.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct UnstrandedFeatureCounts {}

impl Quant for UnstrandedFeatureCounts {
    fn output_only_one_column(&self) -> bool {
        true
    }

    fn weight_read_group(
        &self,
        annotated_reads: &mut [(crate::engine::AnnotatedRead, usize)],
    ) -> Result<()> {
        for (read, _org_index) in annotated_reads.iter_mut() {
            match read {
                crate::engine::AnnotatedRead::Counted(read) => {
                    for (k, v) in read.genes_hit_reverse.drain() {
                        read.genes_hit_correct.insert(k, v); // do not worre about the weight
                    }
                    if read.genes_hit_correct.len() == 1 {
                        //we keep the read, count it as 1.
                    } else {
                        read.genes_hit_correct.clear();
                    }
                }
                _ => {
                    unreachable!();
                }
            }
        }
        Ok(())
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

    fn weight_read_group(
        &self,
        annotated_reads: &mut [(crate::engine::AnnotatedRead, usize)],
    ) -> Result<()> {
        for (read, _org_index) in annotated_reads.iter_mut() {
            match read {
                crate::engine::AnnotatedRead::Counted(info) => {
                    if info.genes_hit_correct.len() == 1 {
                        //we keep the read, count it as 1.
                    } else {
                        *read = crate::engine::AnnotatedRead::FilteredInQuant
                    }
                }
                _ => {
                    unreachable!();
                }
            }
        }
        Ok(())
    }
}
