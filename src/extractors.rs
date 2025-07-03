use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Deserializer};
use serde_valid::Validate;

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
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum UMIExtraction {
    NoUMI(NoUMI),
    #[serde(alias = "regex_name")]
    RegexName(RegexName),

    #[serde(alias = "read_region")]
    ReadRegion(ReadRegion),
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
pub struct RegexName {
    #[serde(deserialize_with = "regex_from_string")]
    regex: regex::Regex,
}

impl UMIExtractor for RegexName {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<String> {
        let name = std::str::from_utf8(read.qname()).ok()?;
        self.regex
            .captures(name)
            .and_then(|caps| caps.get(1).map(|m| m.as_str().to_string()))
    }
}

fn validate_start_stop_rang(start: u16, stop: u16) -> Result<(), serde_valid::validation::Error> {
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
#[validate(custom = |s| validate_start_stop_rang(s.start, s.stop))]
pub struct ReadRegion {
    start: u16,
    stop: u16
}


impl UMIExtractor for ReadRegion {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Option<String> {
        let seq = read.seq().as_bytes()[self.start as usize..self.stop as usize]
            .to_vec();
        Some(String::from_utf8(seq).ok()?)
    }
}
