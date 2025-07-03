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

- do we really need float count values?


- we need a sanity check that no gene interval is disjoint.



- star solo parity:

Zum einen interpretier ich nicht jeden read gleich?
    ff-home:~/upstream/mbf-bam-quantifier/test_cases/starSolo_targeted
    >samtools view actual/output/annotated.bam | grep -P "XQ:[^\t]+" -o | sort | uniq | wc -l
    335
                                                                                                        [ 0s026 | Jul 03 04:27PM ]

     nix develop: /home/user/upstream/mbf-bam-quantifier/test_cases/starSolo_targeted/ 
    ff-home:~/upstream/mbf-bam-quantifier/test_cases/starSolo_targeted
    user>samtools view actual/output/annotated.bam | grep -P "GN:[^\t]+" -o | sort | uniq | wc -l
    345

Das scheinen alles 'schlecht alignierte' Reads zu sein - vermutlich irgendein
filter in STAR den ich übersehe...

das scheinen all

- matrix output isn't stable.
- I'll need to write it to chunk files to enforce an order, I'm afraid.
