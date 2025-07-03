echo "gene_id	count" >../output/counts.tsv
htseq-count \
	-m intersection-nonempty \
	--nonunique all \
	--secondary-alignments score \
	--supplementary-alignments score \
	../input.bam \
	../input.gtf.gz | \
	grep -v "__" \
	>>../output/counts.tsv
