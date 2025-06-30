use std::collections::HashMap;

use super::Quant;
use crate::config::{Input, Output};
use anyhow::bail;
use rust_htslib::bam::{self, Read};
use serde::{Deserialize, Serialize};
use std::io::Write;

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Quantification {}

impl Quantification {
    fn open_output(&self, output: &Output) -> anyhow::Result<std::io::BufWriter<ex::fs::File>> {
        let output_file = format!("{}/counts.tsv", output.directory);
        ex::fs::create_dir_all(&output.directory)?;

        let output_file = ex::fs::File::create(&output_file)?;
        let mut out_buffer = std::io::BufWriter::new(output_file);
        out_buffer
            .write_all(b"reference\tcount\n")
            .expect("Failed to write header to output file");
        Ok(out_buffer)
    }
}

impl Quant for Quantification {
    fn quantify(&mut self, input: &Input, output: &Output) -> anyhow::Result<()> {
        // Implement the quantification logic here
        let bam_with_index = input.get_indexed_bam_reader();
        match bam_with_index {
            Ok(mut bam_with_index) => {
                let index_stats = bam_with_index.index_stats()?;
                let header = bam_with_index.header();
                let reference_names = header.target_names();
                let mut out_buffer = self.open_output(output)?;
                for (tid, length, mapped_count, unmapped_count) in index_stats {
                    if tid == -1 {
                        out_buffer.write(format!("*\t{}\n", unmapped_count).as_bytes())?;
                    } else {
                        let reference_name =
                            std::str::from_utf8(reference_names[tid as usize]).unwrap_or("Unknown");
                        out_buffer
                            .write(format!("{}\t{}\n", reference_name, mapped_count).as_bytes())?;
                    }
                }
            }
            Err(_) => {
                //fall back to counting
                let mut bam = input.get_bam_reader()?;
                let mut counter: HashMap<usize, usize> = HashMap::new();
                let mut unmatched_count = 0;
                let mut out_buffer = self.open_output(output)?;
                let mut read: bam::Record = bam::Record::new();
                while let Some(bam_result) = bam.read(&mut read) {
                    bam_result?;
                    let tid = read.tid();
                    if tid == -1 {
                        unmatched_count += 1;
                        continue;
                    } else {
                        let count = counter.entry(tid as usize).or_insert(0);
                        *count += 1;
                    }
                }
                //now convert tids to contig names
                let header = bam.header();
                let contigs = header.target_names();
                for (reference_sequence_id, contig) in contigs.iter().enumerate() {
                    let count = counter.get(&reference_sequence_id).unwrap_or(&0);
                    let contig = std::str::from_utf8(contig).expect("Reference name wasn't utf-8");
                    out_buffer.write(format!("{}\t{}\n", contig, count).as_bytes())?;
                }
                out_buffer.write(format!("*\t{}\n", unmatched_count).as_bytes())?;
            }
        }
        Ok(())
    }

    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.input.gtf.is_some() {
            bail!("GTF file is not required for contig quantification, do net set it.");
        }

        Ok(())
    }
}
