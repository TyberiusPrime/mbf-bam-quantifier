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


## UMI Grouping   

What makes a UMI a duplicate?

Default is `unique`, which means every UMI is counted (once).

TODO: Umi-tools like options. StarSolo has 1MM_all, 1MM_Directional, 1MM_Directional_UMItools

## Bucket

Conceptually, you only UMI deduplicate reads within one bucket. 

Depending on your sequencing library preparation, you can either 
make use of the read position (`PerPosition`), or count all UMIs 
within one gene (`PerGene` (TODO)).

We also have a mode `PerReference` which is an efficient 'PerGene'
for targeted sequencing aligned to custom genomes. (And that one 
is already implemented).


