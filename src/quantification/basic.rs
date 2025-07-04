use super::{Quant};
use anyhow::Result;



/// count every read that matches. Considers strandedness.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct Basic {
}

impl Quant for Basic {

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
