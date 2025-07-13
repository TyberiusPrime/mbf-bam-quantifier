---
title: Introduction
type: docs
---

# mbf-bam-quantifier

{{< columns >}}
Fast, reliable, flexible region based quantification.

<--->

Count reads just like *you* want.


{{< /columns >}}

Takes a BAM file, a region definition and counts the reads.

Optional unique-molecular-identifier (UMI) based deduplication
and cell barcode based quantification.

## Example 

```toml
[input]
    bam = 'my_aligned_reads.bam'

[input.source]
    mode = 'gtf'
    filename = 'Homo_sapiens.GRCh38.114.chr.gtf.gz'
    feature = exon
    id_attribute = 'gene_id'
    aggr_id_attribute = 'gene_id'

[umi] # where do the umis come from. 
    # If set, trigger UMI deduplication,
    mode = 'tag'
    tag = 'UR'
    grouping = 'unique'
    bucket = 'PerPosition'

[cell_barcodes]
    # if set, trigger single cell output
	extract = {mode = "Tag", tag = "CB"} 
    separator_char = '_'
	max_hamming = 0
    whitelist_files = [
		'barcodes_1.txt',
		'barcodes_2.txt',
		'barcodes_3.txt',
	]

[[filter]] # a list of filters
    mode = 'secondary'
    action = 'remove' # i.e. keep only primary alignments

[strategy]
    direction = 'forward'  # the default


[output]
    directory = 'output'
    write_annotated_bam = false 
```

Run ```mbf-bam-quantifier input.toml```, receive a MTX formated file and statistics, and
an optional annatoted BAM.



