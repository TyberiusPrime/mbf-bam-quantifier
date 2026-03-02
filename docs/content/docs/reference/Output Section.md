---
weight: 4
---

# Output section


```toml
[output]
    directory = "output" # Where do we place the output files?
    write_annotated_bam = false # if set to true, write <directory>/annotated.bam
    mode = Region|SingleCell|StartPositions|Coverage|None # optional, see below.

```

Controls where we place the output.

## Output modes

If you leave mode off, it will default to either 
single cell quantification (if cell_barcodes are present) or region quantification (if no cell_barcodes are present). You can overwrite this if you want to, for example, deduplicate per cell, 
but count per region.

### Region

- a `counts.tsv` file with the counts per region,
  (with columns being <region_id> count_correct count_reverse
- a `counts.tsv.stats.tsv' with some summary statistics

The region id comes from your input definition. For GTF, 
either the aggr_id_attribute or the id_attribute is used,
for references it's 'reference', and for BAM_tags it's the two letter tag.

### SingleCell

(if a `[cell_barcodes]` section is present)

- [Matrix Market Exchange Format](https://math.nist.gov/MatrixMarket/formats.html)
  (that's a matrix.mtx, features.tsv, barcodes.tsv).

### StartPositions 

- a `start_positions.tsv` file with the count of each (corrected) read start postion
   for both strands.
- a `start_positions.tsv.stats.tsv` with some summary statistics

Positions are in genomic coordinates, and 0 based.

Only reads that 'hit' regions are counted. 
Perhaps use bam_references as your region source if you want to count them all.

### Coverage

- a `coverage.tsv` file with the coverage of each detected position.
  Two count columns, so accounting for both strands.

Positions are in genomic coordinates, and 0 based.

Only reads that 'hit' regions are counted. 
Perhaps use bam_references as your region source if you want to count them all.


### None
    Don't create count outputs.

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
- 7 - the read was an approximate UMI duplicate
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


