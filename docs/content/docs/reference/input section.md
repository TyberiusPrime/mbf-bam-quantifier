---
weight: 4
---

# Input section

The input section is mandatory.

```toml
[input]
    bam = "filename.bam"  # What BAM file to quantify
    correct_reads_for_clipping = false # optional see below.
    max_skip_length = 150 # optional see section on correct_reads_for_clipping

# one of
[input.source]
    mode = 'bam_references'
    # uses 0..length, forward for each 'reference' in the BAM file.
    # (see the References filter if you want just some of them)

[input.source]
    # assume the BAM is already annotated with the regions we're supposed to count.
    mode = 'bam_tag'
    tag = 'XX' 

[input.source]
    mode = 'gtf'
    filename` = "something.gtf" # may be gz
    subformat = "GTF|GFF" 
    feature = "exon" # what feature are we counting
    id_attribute = "exon_id" # the attribute to use as the ID of the feature
    aggr_id_attribute = "gene_id" # optional. If set & different from id_attribute, we sum over this.
    duplication_handling = "Collapse|Remove" # what to do with regions that occur multiple times in the GTF
                                             # e.g. exons that occur in multiple transcrpts.
                                             # Defaults to Collapse.
```

## Bam files
BAMs must be sorted and indexed in order to be quantified. 
You can use `samtools sort` and `samtools index` to do this.


## Read positon correction

The argument `correct_reads_for_clipping` controls whether we adjust read positions 
for (soft) clipping, like umi-tools does. Only relevant if you deduplicate PerPosition.

Does not influence the matching of reads to regions (genes).

`max_skip_length` controls the maximum length of a soft clip that you expect to see in your BAM.
We can only 'finish' processing a block of reads that had BAM positions of `x` once we're past `x + max_skip_length`.
You therefore want to keep this lowish to lower memory usage.

Setting this too low (and enabling `correct_reads_for_clipping`) will result in a runtime failure, 
not silently wrong results.


## Sources

### Bam_tag

If you already have a bam tag that defines what region this read belongs to,
you can use the `bam_tag` source.

### GTF

mbf-bam-quantifier can read (GTF and GFF)[https://www.ensembl.org/info/website/upload/gff.html] files (subformat is not not auto-detected at the moment).
Trigger both with `mode = 'GTF'` and the appropriate subformat.

### Duplicate region handling
GTFs that contain transcripts have repeated regions, e.g. exons that occur in multiple transcripts.

By default, we 'collapse' them - identical exons (with identical exon_id) will be counted once,
and will be summed only once into the final gene count.

The other option is to 'remove', which will ignore any reads in these exons
(And report them multiple times in the outut, with 0 reads each!). 
That's what featureCounts does. (I'm not convinced this is actually a useful mode, but 
we want to be able to emulate featureCounts).






