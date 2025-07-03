echo "gene_id	count" >../output/counts.tsv
htseq-count \
	-m intersection-nonempty \
	--nonunique none \
	--secondary-alignments ignore \
	--supplementary-alignments score \
	../input.bam \
	../input.gtf.gz | \
	grep -v "__" \
	>>../output/counts.tsv
