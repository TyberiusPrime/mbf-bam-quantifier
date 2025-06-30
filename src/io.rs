use anyhow::{Context, Result};
use std::io::Read;
use std::path::{Path, PathBuf};

pub fn open_file(filename: impl AsRef<Path>) -> Result<Box<dyn Read + Send>> {
    let fh = ex::fs::File::open(filename.as_ref())
        .context(format!("Could not open file {:?}", filename.as_ref()))?;
    let wrapped = niffler::send::get_reader(Box::new(fh))?;
    Ok(wrapped.0)
}

pub fn open_indexed_bam(
    filename: impl AsRef<Path>,
    index_filename: Option<impl AsRef<Path>>,
) -> Result<rust_htslib::bam::IndexedReader> {
    let mut reader = if let Some(index_filename) = index_filename {
        let filename: &Path = filename.as_ref();
        let filename = filename.to_owned();
        let index_filename_path: PathBuf = index_filename.as_ref().into();
        //htslib wants the same AsRef<Path> twice
        rust_htslib::bam::IndexedReader::from_path_and_index(&filename, &index_filename_path)
            .context(format!("Failed to open indexed BAM file: {:?}", &filename))?
    } else {
        rust_htslib::bam::IndexedReader::from_path(&filename)
            .context(format!("Failed to open BAM file: {:?}", &filename.as_ref()))?
    };
    Ok(reader)
}
