#!/usr/bin/env bash
STAR \
	--runMode genomeGenerate \
	--genomeDir genome.star.index \
	--genomeFastaFiles genome.fasta \
	--sjdbGTFfile genome.gtf \
	--genomeSAindexNbases 7 \
	--sjdbOverhang 100
