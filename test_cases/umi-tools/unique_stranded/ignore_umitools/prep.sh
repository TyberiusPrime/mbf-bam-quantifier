#umi-tools does not quantify by umi and position
# (it's designed for scRnaseq proctolos that first amplify then fragment,
# apparently.
# Suggest method is to use dedup, then count with featureCounts / HTSeq
# but we can ust use our own...

set -eou pipefail
python extract_umi_to_tag.py ../input.bam umitools_input.bam
umi_tools dedup	 -L test.log  --random-seed=123456789 --method=unique --extract-umi-method=tag --stdin=umitools_input.bam --stdout=umitools_dedup.bam --umi-tag=XX --multimapping-detection-method=NH
umi_tools group	 -L test.log  --random-seed=123456789 --method=unique --extract-umi-method=tag --stdin=umitools_input.bam --stdout=umitools_group.bam --umi-tag=XX --multimapping-detection-method=NH --output-bam
samtools index umitools_dedup.bam
samtools index umitools_group.bam
cargo run --release --bin mbf-bam-quantifier quantify-umitools.toml
cp output/counts.tsv ../output/counts.tsv
echo "prepped supposed outputs"

