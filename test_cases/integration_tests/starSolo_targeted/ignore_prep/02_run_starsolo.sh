# swapped R1&R2 is correct.
mkdir -p starsolo_out

	#--soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25, \
	#--soloUMIposition 3_10_3_17 \
STAR \
	--genomeDir genome.star.index \
	--genomeLoad NoSharedMemory \
	--soloType CB_UMI_Complex \
	--soloCBwhitelist ../input_BD_CLS1.txt ../input_BD_CLS2.txt ../input_BD_CLS3.txt \
	--readFilesIn "input_R2.fasta" "input_R1.fasta" \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--soloCBposition "0_0_0_8" "0_21_0_29" "0_43_0_51" \
	--soloUMIposition "0_52_0_60" \
	--soloCBmatchWLtype Exact \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMultimapScoreRange 0 \
	--seedSearchStartLmax 50 \
	--clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
	-- outFileNamePrefix starsolo_out/

