use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Deserializer};

pub fn regex_from_string<'de, D>(deserializer: D) -> core::result::Result<regex::Regex, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    let re = regex::Regex::new(&s)
        .map_err(|e| serde::de::Error::custom(format!("Invalid regex: {e}")))?;
    Ok(re)
}
#[enum_dispatch(UMIExtraction)]
pub trait UMIExtractor {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<String>;
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
#[serde(tag = "type")]
#[enum_dispatch]
pub enum UMIExtraction {
    NoUMI(NoUMI),
    #[serde(alias = "regex_name")]
    RegexName(UMIExtractionRegexName),
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
    fn extract(&self, _read: &rust_htslib::bam::record::Record) -> Option<String> {
        None // No UMI extraction, return None
    }
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct UMIExtractionRegexName {
    #[serde(deserialize_with = "regex_from_string")]
    regex: regex::Regex,
}

impl UMIExtractor for UMIExtractionRegexName {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<String> {
        let name = std::str::from_utf8(read.qname()).ok()?;
        self.regex
            .captures(name)
            .and_then(|caps| caps.get(1).map(|m| m.as_str().to_string()))
    }
}
