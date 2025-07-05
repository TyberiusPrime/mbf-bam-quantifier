- smart exon tag counting like we do...

- more umitools comparisons
- more htseq comparisons.

- star solo comparison tests

- coverage quantification?

- faster interval format instead of gtf parsing?

- can we split our intervals bgfz borders? 


- options to skip filtered / duplicates on write_output_bam, 
  for umi-tools feature parity?

- htseq also has a scrnaseq counter mode...


- what happens whith multi mappers that map multiple times in one read
    - my old code didn't count the multiple times. think we should reimplement that...
    - perhaps by a filter on the AnnotatedReads...

- gene level quantification for the non -position specific scRNAseqs..

- write some documentation

- add sam header for this program

- do we really need float count values?

- we need a sanity check that no gene interval is disjoint / occurs on multiple references.

- star solo parity, for multi-mappers

- add test case for reference filter

- fancy umi clustering algorithms like umi-tools, (and like starsolo implements them)
- add individual read consideration test cases (the compare-to-others approach is difficult to debug when it fails).
- add star solo spliced reads test case (and dataset)
- figure out why starsolo is padding / truncated umis to 10bp (what?!)
- memory limitation workarounds.
- do the umi dedup before assign the reads to the genes - that should give us a bit of a speed boost / same some ram/


- I think I can use less memory and save on sorting, if I 
  change the annotate-reads->sort-reads->count to a
  'count reads until next position (that takes a &mut HashMap<(position) ,annotated-reads-forward/reverse), 
  then when that returns, processes all positions which no position correction for S 
  can possibly occur, and then discards those eagerly. add in a final flush, and that should work,
  but limit our memory usage from <all-reads-within-chunk> to <max-no-of-reads-in-a-max-skip-sized-area>,
  the later being always smaller

