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
  (problem if doing umi-dedup per gene and not reference=gene)

- test case for umi-dedup-per-reference

- umi dedup per gene (not reference!)
- umi dedup per (start, stop)?
- test case for umi-n-filter


- starsolo test case that shows you need                
    {"mode": "NInUMI", "action": "remove"},
    {"mode": "PolyAUMI", "action": "remove"}, #what are you doing STAR?


- I'm not happy with the configuration. singlecell output should'nt be in dedup,
dedup + umi should be one...

- compressed output (gzip/zstd)

- zstd compressed gtfs
- is the id_attribute aggr_id_attribute split sensible? add explanation 
  I think there's a bug lurking here if multi-region and id=attr_id, we're not splitting at the right position?
  And id != attr_id actually fails (no test case...)
  Add test case, rework setting


- the genes hit storage doesn't belong on the reads,
  but on the umi-dedup bucket! 

+- checkout https://github.com/Daniel-Liu-c0deb0t/UMICollapse
 (doesn't seem to have test cases)


umi dedup methods
- [x] unique
- [ ] percentile
- [ ] cluster, 
- [ ] adjacency 
- [ ] directional,
- [ ] starsolo directional,
- [ ] BD DBEC



- when deduping per gene/reference, we need to make sure correct_read_pos is false.
Otherwise we might process the fake '0' position somewhere after max_split_size again.
(our paranoia code should catch that though, so we'll know when it happens)
