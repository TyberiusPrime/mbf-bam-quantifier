# mbf-bam-quantifier

A most flexible, structured bam quantifier.

It filters reads, deduplicates on UMIs, counts in regions or genes,
for bulk or single cell sequencing, stranded or unnstranded.

You get to choose 
 - where your regions come from
 - what reads to consider,
 - where barcodes and UMIs are being read from,
 - what's considered a UMI duplicate, b
 - how and if barcodes are matched to a whit list and 
 - whether reads must match the region's strand.


It's fast (~ 100 million reads a minute, give or take)
and has extensive test coverage.





