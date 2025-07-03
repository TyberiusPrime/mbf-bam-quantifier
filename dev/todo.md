- smart exon tag counting like we do...

- umitools comparisons

- star solo comparison tests

- coverage quantification?

- faster interval format instead of gtf parsing?

- can we split our intervals bgfz borders? 


- options to skip filtered / duplicates on write_output_bam, 
  for umi-tools feature parity?

- compare to htseq
- htseq also has a scrnaseq counter mode...

 - paired end handling
   at least a filter read2 mode.

 -rename filter mode to type, or make it consistnent across all config

- what happens whith multi mappers that map multiple times in one read
    - my old code didn't count them. think we should reimplement that...
    - perhaps by a filter on the AnnotatedReads...
