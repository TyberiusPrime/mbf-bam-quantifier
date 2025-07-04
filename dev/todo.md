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

Moeglich das Star will das die reads vollstaendig im transkript liegen?
intersection-strict doesn't change a thing.
Why would it, we're looking at gene-is-contig bams.

I think there are two cases.
NB552003:178:HWGJ7BGXL:1:11101:25683:1808 mappes 4 times, and starsolo counts it on one place,
and doesn't count the others. (soloMultiMappers Unique, I suppose).
Of course, having the *unsorted* reads at hand, it's in a good position to do that.
While I'm not...
But that's not what 'unique' is supposed to do, I think. It should not count the read 
*at all*. Uniform also doesn't change it to be multi. 


NB552003:178:HWGJ7BGXL:1:11101:13470:1812	16	Cxcl5|ENSMUST00000031318.5|Reference_end	130	255	49S25M2S	*	0	0	ATCAGTTTTTTATCTTTTGGAATCCCTGCTTACAAGTCTTGAAGAAAATCATGTTTGATTGTGTAGTGTTGTGACC	///////E/E<////EE/////E///<//E//////E/EE////E///////E6E///E/EEEEE/E/E/AAA///	NH:i:1	HI:i:1	nM:i:6	AS:i:12	CR:Z:GCGATTACA_ACGATGAAT_ACAGTAAAC	UR:Z:TTAAACAAT	GX:Z:-	GN:Z:-	sS:Z:GCGATTACAACTGGCCTGCGAACGATGAATGGTAGCGGTGACAACAGTAAACTTAAACAATTTTTTTTTCTTTTTT	sQ:Z:AAAAAEEEEEEAEEEEEEEEE/E/EE/EEEEEEEEEEEEEEEEE/EE/EE/EEEE6E<EEEEEEE6A<E/<EAAAE	sM:i:0	CB:Z:GCGATTACA_ACGATGAAT_ACAGTAAAC	UB:Z:ATTAAACAAT


Maybe the whole multi-gene thing is a red herring

why is it padding the UMI?
Always forcing it to be 10 bp long???
actually, it's also truncated longer umis to 10b...
that's gotta be a bug?

Even without mult mappers, I assign genes where starsolo doesn.t
NB552003:178:HWGJ7BGXL:1:11101:2082:1632	16	Lipg|ENSMUST00000066532.4|Reference_end	246	255	31S19M26S	*	0	0	TTCATTAATTTTGCTTGTTTGTACTTAGTTTTTTTTCATTTAAGATTAACATTCATTATTATGCTGTGGCCTTTTC	EEE//E//E/EE6/A//EE//AE/EEE/EEEAE/EA/AA//A/AA/A/E/E6E//EEE6//E/EE///AE//AA//	NH:i:1	HI:i:1	nM:i:2	AS:i:14	CR:Z:TATGTGGCA_TTGCGTACA_TTCAGCTCA	UR:Z:GTACGGTTT	GX:Z:-	GN:Z:-	sS:Z:TATGTGGCAACTGGCCTGCGATTGCGTACAGGTAGCGGTGACATTCAGCTCAGTACGGTTTTTTTTTTTTTTTTTT	sQ:Z:/AAAAEEEAEEEEEEEEEEEEEEEEE/AEEEEEEEEEEEEEEEEEEEEEAEEE/EE/6/EEEEAEAEAEEEEAE<E	sM:i:0	UB:Z:AGTACGGTTT	XQ:Z:Lipg|ENSMUST00000066532.4|Reference_end=1.00	XR:Z:	XP:i:246	CB:Z:TATGTGGCA_TTGCGTACA_TTCAGCTCA


maybe it's strand mismatch... but I also have reads 
that ain't reversed

it only counts them iff they agree with the transcript's splicing (from the GTF)









das scheinen all

- matrix output isn't stable.
- I'll need to write it to chunk files to enforce an order, I'm afraid.
