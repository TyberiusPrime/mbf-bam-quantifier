#count dedupes per gene, not per read position group... 
#umi_tools count -L test.log  --random-seed=123456789 --method=directional --gene-tag=XF --skip-tags-regex="^[__|Unassigned]" --extract-umi-method=umis --stdin=input_chr19_gene_tags.bam
#umi_tools count -L test.log  --random-seed=123456789 --method=unique --gene-tag=XF --skip-tags-regex="^[__|Unassigned]" --extract-umi-method=umis --stdin=input_chr19_gene_tags.bam
#umi_tools group	 -L test.log  --random-seed=123456789 --method=unique --extract-umi-method=umis --stdin=input_chr19_gene_tags.bam --stdout=ignore_umitools_group.bam --output-bam

# this is where you need to coun the reads in
umi_tools dedup	 -L test.log  --random-seed=123456789 --method=unique --extract-umi-method=umis --stdin=input_chr19_gene_tags.bam --stdout=ignore_umitools_dedup.bam
samtools view ignore_umitools_dedup.bam | grep -P ENSG[0-9]\+ -o | sort | uniq -c
