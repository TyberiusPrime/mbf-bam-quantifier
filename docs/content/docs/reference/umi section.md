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
    external_umi_thresholder_command = ["umi_threshold.py"] # optional. see below.
```

## Bucket

Conceptually, you only UMI deduplicate reads within one bucket.

Depending on your sequencing library preparation, you can either
make use of the read position (`PerPosition`), or count all UMIs
within one region (`PerRegion`).

Note that with PerRegion, your regions must not overlap.
(It's fine to have overlapping sub-regions within one regions in your GTF,
but the regions we ultimately quantify UMIs in must not overlap).

We also have a mode `PerReference` which is an efficient 'PerGene'
for targeted sequencing aligned to custom genomes.

If you're doing single cell (i.e. you have cell barcodes),
buckets are further split by (corrected) cell barcode.

## UMI Grouping

What makes a UMI a duplicate?

### unique

Default is `unique`, which means every UMI is counted (once).

### percentile

UMIS that have a count below 1% of the median of the UMis within the bucket are considered duplicates.

(umi-tools: [`percentile`](https://umi-tools.readthedocs.io/en/latest/the_methods.html))

### Cluster

All UMIs that are within `max_distance = ...` (Hamming) of each other are considered duplicates.
(The chosen read for the annotated.bam is one from the UMI with the highest count).

(umi-tools: [`cluster`](https://umi-tools.readthedocs.io/en/latest/the_methods.html), STARSolo: `1MM_all`)

### Directional

Form networks with edges defined based on hamming distance threshold
(`max_distance = ...`) and `node A counts >= (2 * node B counts) - 1`.

Each connected component is considered a UMI group.

(umi-tools: [`1MM_Directional`](https://umi-tools.readthedocs.io/en/latest/the_methods.html),
STARSolo: `1MM_Directional_UMItools`)

BD Rhapsody calls this method RSEC.

### Directional_STARSolo

Form networks with edges defined based on hamming distance threshold
(`max_distance = ...`) and `node A counts >= (2 * node B counts) - 0`.

(Note the offset difference to 'Directional').

Each connected component is considered a UMI group.

Described by STARSolo as a "same as 1MM_Directional_UMItools [see above], but
with more stringent criteria for duplicate UMIs".

I believe, this essentialy just disconnects UMIs that have a count of 1 from
other UMIs in the cluster that also have a count of 1.

(StarSolo: `1MM_Directional`)

### Todo

StarSolo 1MM_all, umi-tools adjancent (same thing, mabye?)

## external_umi_thresholder_command

Some sequencing methodologies suggest running a
thresholding on UMIs per biomolecule (e.g. DBEC in BD Rhapsody(TM)).

Since these can be arbitrarily complex, we outsource these to external
scripts, which you can specify with `external_umi_thresholder_command` (as a
list of command line arguments).

The receiving command will receive a comma-separated-list of umi-counts
per biomolecule on stdin.

It is supposed to return a single (integer, > 0) number on stdout:
The threshold below which we'll remove the UMIs from counting.

Returning anything else on stdout will result in an error
(and an aborted quantification). Output on stderr is ignored (but shown in case of a
non-zero return code / a non-parsable result).
