---
weight: 4
---

# Output section


```toml
[output]
    directory = "output" # Where do we place the output files?
    write_annotated_bam = false # if set to true, write <directory>/annotated.bam

```

Controls where we place the output.

## Output for non single cell quantification

(if no `[cell_barcodes]` section is present)


- a `counts.tsv` file with the counts per region,
  (with columns being <region_id> count_correct count_reverse
- a `counts.tsv.stats.tsv' with some summary statistics

The region id comes from your input definition. For GTF, 
either the aggr_id_attribute or the id_attribute is used,
for references it's 'reference', and for BAM_tags it's the two letter tag.

## Output for single cell quantification

(if a `[cell_barcodes]` section is present)

- [Matrix Market Exchange Format](https://math.nist.gov/MatrixMarket/formats.html)
  (that's a matrix.mtx, features.tsv, barcodes.tsv).

## annotated.bam

If requested, we output the decisions on each read as `<directory>/annotated.bam`.

We add the following tags (and remove their old values if they were set in the BAM file):

Note that depending on where exactly reads are filtered or detected as duplicates,
some of the tags may not be set.

### XF:i - filter decision

- 1 - the read was removed by a filter
- 3 - the read was detected as an UMI duplicate
- 4 - the read's cell barcode  was not in the whitelist
- 5 - the read had no barcode
- 6 - the read had no UMI
- counted reads do not get an XF tag.

### XQ - correct hits

Genes (regions) hit in the correct orientation, comma separated 

### XR - reverse hits
Genes (regions) hit in the incorrect orientation, comma separated 

### XP - corrected position

The corrected position of the read. See TODO.
Only present if input.correct_reads_for_clipping is set to true.

### CB corrected cell barcode
The (corrected) cell barcode.

### CR uncorrected cell barcode
If a read's barcode was not in the white list, the uncorrected barcode is stored here.
(otherwise this tag is not set)

### UMI CR / UR
TODO (not yet implemented)


