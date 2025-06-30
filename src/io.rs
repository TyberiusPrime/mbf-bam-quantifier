use anyhow::{Context, Result};
use rust_htslib::bam;
use std::io::Read;
use std::path::Path;

pub fn open_file(filename: impl AsRef<Path>) -> Result<Box<dyn Read + Send>> {
    let fh = ex::fs::File::open(filename.as_ref())
        .context(format!("Could not open file {:?}", filename.as_ref()))?;
    let wrapped = niffler::send::get_reader(Box::new(fh))?;
    Ok(wrapped.0)
}

pub fn open_indexed_bam(
    bam_filename: impl AsRef<Path>,
    bai_filename: Option<impl AsRef<Path>>,
) -> Result<bam::IndexedReader> {
    let bam = match bai_filename {
        Some(ifn) => bam::IndexedReader::from_path_and_index(bam_filename.as_ref(), ifn.as_ref()),
        _ => bam::IndexedReader::from_path(bam_filename.as_ref()),
    };
    bam.context("Could not read bam: {}")
}
