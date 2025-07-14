---
title: Introduction
type: docs
---

# mbf-bam-quantifier

{{< columns >}}
Fast, reliable and flexible region based quantification.

<--->

Count reads just like *you* want.


{{< /columns >}}

Take a BAM file, a region definition and counts the reads.

Optional unique-molecular-identifier (UMI) based de-duplication
and cell barcode based quantification.

## Example 

```toml
[input]
    bam = 'my_aligned_reads.bam'

[input.source]
    mode = 'gtf'
    filename = 'Homo_sapiens.GRCh38.114.chr.gtf.gz' # e.g. from ensembl
    feature = exon
    id_attribute = 'gene_id'

[umi] 
    # If set, trigger UMI deduplication,
    # How do we group umis
    grouping = 'unique'
    bucket = 'PerPosition'

    # where do they come from?
    extract = {mode = 'tag', tag = 'UR'}

[cell_barcodes]
    # if set, trigger single cell output
	extract = {mode = "Tag", tag = "CB"} 
	max_hamming = 0

    separator_char = '_'
    whitelist_files = [ 'barcodes_1.txt', 'barcodes_2.txt', ]

[[filter]] # zero or more filters to the reads.
    mode = 'secondary'
    action = 'remove' # i.e. keep only primary alignments

[strategy]
    direction = 'forward'  # only reads matching the gtf direction.


[output]
    directory = 'output'
    write_annotated_bam = false 
```

Run ```mbf-bam-quantifier input.toml```, receive a MTX formatted file and statistics, and
an optional annatoted BAM.



