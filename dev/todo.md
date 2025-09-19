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
- [x] percentile
- [ ] cluster, 
- [ ] adjacency 
- [ ] directional,
- [ ] starsolo directional,
- [ ] BD DBEC



- when deduping per gene/reference, we need to make sure correct_read_pos is false.
Otherwise we might process the fake '0' position somewhere after max_split_size again.
(our paranoia code should catch that though, so we'll know when it happens)

- failure test case for gtf that looks like this
1	mbf	tss	16437	18436	0	.	0	gene_stable_id "ENSG00000278267"
1	mbf	tss	28571	30570	0	.	0	gene_stable_id "ENSG00000227232"
1	mbf	tss	35082	37081	0	.	0	gene_stable_id "ENSG00000237613"
1	mbf	tss	51473	53472	0	.	0	gene_stable_id "ENSG00000268020"
1	mbf	tss	56598	58597	0	.	0	gene_stable_id "ENSG00000290826"
1	mbf	tss	61949	63948	0	.	0	gene_stable_id "ENSG00000240361"
1	mbf	tss	64419	66418	0	.	0	gene_stable_id "ENSG00000186092"
1	mbf	tss	90106	92105	0	.	0	gene_stable_id "ENSG00000239945"
1	mbf	tss	130025	132024	0	.	0	gene_stable_id "ENSG00000233750"
1	mbf	tss	132724	134723	0	.	0	gene_stable_id "ENSG00000238009"

(missing final ;)

- what's the unstranded / inverse stranded story.
We have that in the strategy.

- document the id_attribute, aggr_id_attribute functionality.
id_attribute -> used in region dedup / featureCounts parity.
aggr_id_attribute -> used to define what we add up, and where the splits can occur.


- profile nested_intervals vs rust bio  intervaltrees
