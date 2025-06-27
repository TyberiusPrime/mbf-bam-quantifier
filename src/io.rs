use anyhow::{Context, Result};
use std::path::Path;
use std::io::Read;


pub fn open_file(filename: impl AsRef<Path>) -> Result<Box<dyn Read + Send>> {
    let fh = ex::fs::File::open(filename.as_ref())
        .context(format!("Could not open file {:?}", filename.as_ref()))?;
    let wrapped = niffler::send::get_reader(Box::new(fh))?;
    Ok(wrapped.0)
}
