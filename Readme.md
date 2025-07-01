# featureCounts

mbf-bam-quantifier can achieve the same read counts as featureCounts.

To do so, 
 - set the quant mode to 'unstranded_featurecounts' or 'stranded_featurecounts',
 - set input.source.duplicate_handling = 'rename'
 - set a filter against multimapping reads.

(you can see the testcases in test_cases/integration_tests/feature_counts_exon_level_stranded/
for an illustration)

This works then the same way as featureCounts, which reports
repeated exons in the source file multiple times, but sets all their reads to 0.

**That's likely not what you want. Or would have expected**




