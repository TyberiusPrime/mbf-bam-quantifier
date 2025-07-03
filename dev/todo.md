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
    - my old code didn't count them. think we should reimplement that...
    - perhaps by a filter on the AnnotatedReads...

- gene level quantification for the non -position specific scRNAseqs..

- write some documentation


- add sam header for this program


 - we need to be able to pus errors from the worker threads
