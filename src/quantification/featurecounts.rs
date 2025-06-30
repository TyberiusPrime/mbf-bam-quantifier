mod chunked_genome;
mod counters;


use super::{Quant, build_trees_from_gtf};
use crate::config::{Input, Output};
use crate::gtf::GTFEntrys;
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Debug, Clone, Serialize)]
pub struct Quantification {
    feature: String,
    id_attribute: String,
}

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
        let bam_with_index = input.get_indexed_bam_reader()?;

        let mut gtf_entrys = if self.feature == "gene" {
            input.read_gtf(&[self.feature.as_str()])?
        } else {
            input.read_gtf(&["gene", &self.feature])?
        };

        let gene_trees =
            build_trees_from_gtf(&self.id_attribute, gtf_entrys.remove("gene").context("No GTF entries found for `gene`s. Necessary to know where to split the genome productivly.")?)
                .context("Failed to build trees from GTF")?;
        let feature_trees = if self.feature == "gene" {
            gene_trees.clone()
        } else {
            build_trees_from_gtf(
                &self.id_attribute,
                gtf_entrys
                    .remove(&self.feature)
                    .context("No GTF entries found for the specified feature")?,
            )?
        };

        let counts = counters::count_reads_stranded(
            &input.bam,
            None,
            gene_trees,
            feature_trees,
            true,
            None,
        );


        Ok(())
    }

    fn check(&self, config: &crate::config::Config) -> anyhow::Result<()> {
        if config.input.gtf.is_none() {
            bail!("GTF file is required for featurecount quantification");
        }

        Ok(())
    }
}
