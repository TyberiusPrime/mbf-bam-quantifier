use anyhow::{bail, Context, Result};
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
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Result<Option<Vec<u8>>>;
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum UMIExtraction {
    #[serde(alias = "regex_name")]
    RegexName(RegexName),

    #[serde(alias = "search_in_name")]
    SearchInName(SearchInName),

    #[serde(alias = "read_region")]
    ReadRegion(ReadRegion),
    #[serde(alias = "Tag")]
    Tag(Tag),
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct RegexName {
    #[serde(deserialize_with = "u8_regex_from_string")]
    regex: regex::bytes::Regex,
}

impl UMIExtractor for RegexName {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Result<Option<Vec<u8>>> {
        let name = read.qname();
        self.regex
            .captures(name)
            .map(|caps| caps.get(1).map(|m| name[m.start()..m.end()].to_vec()))
            .with_context(|| {
                format!(
                    "Failed to extract from read name via regex: {}",
                    std::str::from_utf8(read.qname()).unwrap_or("<invalid UTF-8>")
                )
            })
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
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Result<Option<Vec<u8>>> {
        if self.stop > read.seq_len() as u16 {
            bail!(
                "Read region stop ({}) exceeds read length ({})",
                self.stop,
                read.seq_len()
            );
        }
        let seq = read.seq().as_bytes()[self.start as usize..self.stop as usize].to_vec();
        Ok(Some(seq))
    }
}
#[derive(serde::Deserialize, Debug, Clone, Validate)]
#[serde(deny_unknown_fields)]
pub struct Tag {
    #[serde(deserialize_with = "crate::config::deser_tag")]
    pub tag: [u8; 2],
}

impl UMIExtractor for Tag {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Result<Option<Vec<u8>>> {
        let tag = read.aux(&self.tag)?;
        match tag {
            rust_htslib::bam::record::Aux::String(v) => Ok(Some(v.as_bytes().to_vec())),
            _ => bail!("Expected tag to be a string, found {:?}", tag),
        }
    }
}

#[derive(serde::Deserialize, Debug, Clone)]
#[serde(deny_unknown_fields)]
pub struct SearchInName {
    #[serde(deserialize_with = "crate::config::u8_from_string")]
    search: Vec<u8>,
    skip: usize,
    len: usize,
}

impl UMIExtractor for SearchInName {
    fn extract(&self, read: &rust_htslib::bam::record::Record) -> Result<Option<Vec<u8>>> {
        let name = read.qname();
        //where is search in name?
        let offset = twoway::find_bytes(name, &self.search);
        match offset {
            Some(offset) => {
                let start = offset + self.search.len() + self.skip;
                let stop = start + self.len;
                if stop > name.len() {
                    Ok(None)
                } else {
                    Ok(Some(name[start..stop].to_vec()))
                }
            }
            None => Ok(None),
        }
    }
}
