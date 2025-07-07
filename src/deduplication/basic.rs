use super::Dedup;

/// count every read that matches. Considers strandedness.
#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
#[serde(deny_unknown_fields)]
pub struct Basic {}

impl Dedup for Basic {}
