---
weight: 100
---

# UMI Section

```toml
[umi]
    strategy.kind = 'cluster' # When to consider UMIs duplicates. see Below
    strategy.max_distance = 1
    bucket = 'PerPosition', # What reads to consider when grouping UMIs. see below
    extract = ... # An Extractor. See Extractors in the navigation
```


## Bucket

Conceptually, you only UMI deduplicate reads within one bucket. 

Depending on your sequencing library preparation, you can either 
make use of the read position (`PerPosition`), or count all UMIs 
within one gene (`PerGene` (TODO)).

We also have a mode `PerReference` which is an efficient 'PerGene'
for targeted sequencing aligned to custom genomes. (And that one 
is already implemented).

If you're doing single cell (i.e. you have cell barcodes), 
buckets are further split by (corrected) cell barcode.


## UMI Grouping   

What makes a UMI a duplicate?

Default is `unique`, which means every UMI is counted (once).

The further options are 

## percentile

UMIS that have a count below 1% of the median of the UMis within the bucket are considered duplicates.

(umi-tools: `percentile`)

## Cluster

All UMIs that are within <max-hamming-distance> of each other are considered duplicates.
(The chosen read for the annotated.bam is one from the UMI with the highest count).

(umi-tools: `cluster`, STARSolo: `1MM_all`)

# 1MM_Directional

Form networks with edges defined based on hamming distance threshold and `node A counts >= (2 * node B counts) - 1`.

Each connected component is considered a UMI group.

TODO: implement

(umi-tools: `1MM_Directional`, STARSolo: `1MM_Directional_UMItools`)

# 1MM_Directional STAR

TODO: Figure out what star is doing...


    

TODO: Umi-tools like options. StarSolo has 1MM_all, 1MM_Directional, 1MM_Directional_UMItools
