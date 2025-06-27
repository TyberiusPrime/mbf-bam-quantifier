use std::collections::HashMap;

use super::Quant;
use crate::config::{Input, Output};
use anyhow::{bail, Result};
use noodles_csi::binning_index::ReferenceSequence;
use noodles_csi::BinningIndex;
use serde::{Deserialize, Serialize};
use std::io::Write;

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Quantification {}

impl Quantification {
    fn open_output(&self, output: &Output) -> anyhow::Result<std::io::BufWriter<ex::fs::File>> {
        let output_file = format!("{}/counts.tsv", output.directory);
        ex::fs::create_dir_all(&output.directory)?;
        //what's in .
        for file in std::fs::read_dir(&output.directory)? {
            let file = file?;
            dbg!(file.file_name());
        }

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
                let mut out_buffer = self.open_output(output)?;
                let bai = input.get_bam_index()?;
                dbg!(bai.unplaced_unmapped_record_count());

                let header = bam_with_index.read_header()?;
                let reference_sequence_names = header
                    .reference_sequences()
                    .iter()
                    .map(|(name, _)| std::str::from_utf8(name).unwrap().to_string())
                    .collect::<Vec<_>>();
                //let bai = bam_with_index.index();
                dbg!(&reference_sequence_names);

                for (ii, reference_seq_info) in bai.reference_sequences().iter().enumerate() {
                    let matched_count = reference_seq_info
                        .metadata()
                        .map(|metadata| metadata.mapped_record_count())
                        .unwrap_or(0);
                    out_buffer.write(
                        format!(
                            "{}\t{}\n",
                            reference_sequence_names[ii],
                            matched_count
                        )
                        .as_bytes(),
                    )?;
                }

                let unmatched_count = bai.unplaced_unmapped_record_count().unwrap_or(0);
                out_buffer.write(format!("*\t{}\n", unmatched_count).as_bytes())?;
            }
            Err(_) => {
                //fall back to counting
                let mut counter: HashMap<usize, usize> = HashMap::new();
                let mut unmatched_count = 0;
                let mut reader = input.get_bam_reader()?;
                let mut out_buffer = self.open_output(output)?;
                let header = reader.read_header()?;

                for result in reader.records() {
                    let record = result?;
                    let tid = record.reference_sequence_id();
                    match tid {
                        None => {
                            unmatched_count += 1;
                            continue;
                        }
                        Some(Ok(tid)) => {
                            let count = counter.entry(tid).or_insert(0);
                            *count += 1;
                        }
                        Some(Err(e)) => {
                            bail!("Error getting reference sequence ID: {}", e);
                        }
                    }
                }
                //now convert tids to contig names
                let contigs: Result<Vec<String>> = header
                    .reference_sequences()
                    .iter()
                    .map(|(seq, _)| Ok(std::str::from_utf8(seq)?.to_string()))
                    .collect();
                let contigs = contigs?;

                for (reference_sequence_id, contig) in contigs.iter().enumerate() {
                    let count = counter.get(&reference_sequence_id).unwrap_or(&0);
                    out_buffer
                        .write(format!("{}\t{}\n", contig, count).as_bytes())?;
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
