use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Deserializer};
use serde_valid::Validate;

pub fn u8_regex_from_string<'de, D>(
    deserializer: D,
) -> core::result::Result<regex::bytes::Regex, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    let re = regex::bytes::Regex::new(&s)
        .map_err(|e| serde::de::Error::custom(format!("Invalid regex: {e}")))?;
    Ok(re)
}

#[enum_dispatch(UMIExtraction)]
pub trait UMIExtractor {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<Vec<u8>>;
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum UMIExtraction {
    NoUMI(NoUMI),
    #[serde(alias = "regex_name")]
    RegexName(RegexName),

    #[serde(alias = "read_region")]
    ReadRegion(ReadRegion),
    #[serde(alias = "Tag")]
    Tag(Tag),
}

impl Default for UMIExtraction {
    fn default() -> Self {
        UMIExtraction::NoUMI(NoUMI {})
    }
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct NoUMI {}

impl UMIExtractor for NoUMI {
    fn extract(&self, _read: &rust_htslib::bam::record::Record) -> Option<Vec<u8>> {
        None // No UMI extraction, return None
    }
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct RegexName {
    #[serde(deserialize_with = "u8_regex_from_string")]
    regex: regex::bytes::Regex,
}

impl UMIExtractor for RegexName {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<Vec<u8>> {
        let name = read.qname();
        self.regex
            .captures(name)
            .and_then(|caps| caps.get(1).map(|m| name[m.start()..m.end()].to_vec()))
    }
}

fn validate_start_stop_range(start: u16, stop: u16) -> Result<(), serde_valid::validation::Error> {
    if start >= stop {
        Err(serde_valid::validation::Error::Custom(
            "start must be less than stop".to_string(),
        ))
    } else {
        Ok(())
    }
}

#[derive(serde::Deserialize, Debug, Clone, Validate)]
#[serde(deny_unknown_fields)]
#[validate(custom = |s| validate_start_stop_range(s.start, s.stop))]
pub struct ReadRegion {
    start: u16,
    stop: u16,
}

impl UMIExtractor for ReadRegion {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<Vec<u8>> {
        let seq = read.seq().as_bytes()[self.start as usize..self.stop as usize].to_vec();
        Some(seq)
    }
}
#[derive(serde::Deserialize, Debug, Clone, Validate)]
#[serde(deny_unknown_fields)]
pub struct Tag {
    #[serde(deserialize_with = "crate::config::deser_tag")]
    pub tag: [u8; 2],
}

impl UMIExtractor for Tag {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<Vec<u8>> {
        let tag = read.aux(&self.tag).ok()?;
        match tag {
            rust_htslib::bam::record::Aux::String(v) => Some(v.as_bytes().to_vec()),
            _ => None,
        }
    }
}
