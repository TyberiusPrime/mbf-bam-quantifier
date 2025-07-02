use enum_dispatch::enum_dispatch;

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize, PartialEq,Eq)]
enum KeepOrRemove {
    #[serde(alias = "keep")]
    Keep,
    #[serde(alias = "remove")]
    Remove,
}

#[enum_dispatch(Filter)]
pub trait ReadFilter: Send + Sync {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool;
}

#[derive(serde::Deserialize, Debug, Clone, strum_macros::Display, serde::Serialize)]
#[serde(tag = "mode")]
#[enum_dispatch]
pub enum Filter {
    #[serde(alias = "multimapper")]
    MultiMapper(MultiMapper),
    #[serde(alias = "non_primary")]
    NonPrimary(NonPrimary),
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct MultiMapper {
    action: KeepOrRemove,
}

impl ReadFilter for MultiMapper {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        use rust_htslib::bam::record::Aux;
        let hit = match read.aux(b"NH") {
            Ok(Aux::I8(value)) => value > 1,
            Ok(Aux::I16(value)) => value > 1,
            Ok(Aux::I32(value)) => value > 1,
            Ok(Aux::U8(value)) => value > 1,
            Ok(Aux::U16(value)) => value > 1,
            Ok(Aux::U32(value)) => value > 1,

            _ => {
                panic!(
                    "read had no NH tag (or wasn't an integer). Can't remove multi mappers. Aborting"
                )
            }
        };
        match self.action {
            KeepOrRemove::Keep => !hit,
            KeepOrRemove::Remove => hit,
        }
    }
}

#[derive(serde::Deserialize, Debug, Clone, serde::Serialize)]
pub struct NonPrimary {
    action: KeepOrRemove,
}

impl ReadFilter for NonPrimary {
    fn remove_read(&self, read: &rust_htslib::bam::record::Record) -> bool {
        if read.is_secondary() {
            return self.action == KeepOrRemove::Remove;
        }
        false
    }
}
