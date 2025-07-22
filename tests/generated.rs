
mod test_runner;
use test_runner::run_test;

#[test]
fn test_case_basic_basic_stranded_gene_is_contig() {
    run_test(std::path::Path::new("test_cases/basic/basic_stranded_gene_is_contig"));
}

#[test]
fn test_case_basic_basic_stranded_gene_is_contig_reverse_direction() {
    run_test(std::path::Path::new("test_cases/basic/basic_stranded_gene_is_contig_reverse_direction"));
}

#[test]
fn test_case_basic_basic_unstranded_gene_is_contig() {
    run_test(std::path::Path::new("test_cases/basic/basic_unstranded_gene_is_contig"));
}

#[test]
fn test_case_basic_basic_unstranded_gene_is_contig_incl_multimapper() {
    run_test(std::path::Path::new("test_cases/basic/basic_unstranded_gene_is_contig_incl_multimapper"));
}

#[test]
fn test_case_featurecounts_feature_counts_exon_level_stranded() {
    run_test(std::path::Path::new("test_cases/featurecounts/feature_counts_exon_level_stranded"));
}

#[test]
fn test_case_featurecounts_feature_counts_exon_level_stranded_reversed() {
    run_test(std::path::Path::new("test_cases/featurecounts/feature_counts_exon_level_stranded_reversed"));
}

#[test]
fn test_case_featurecounts_feature_counts_exon_level_unstranded() {
    run_test(std::path::Path::new("test_cases/featurecounts/feature_counts_exon_level_unstranded"));
}

#[test]
fn test_case_featurecounts_feature_counts_gene_level_stranded() {
    run_test(std::path::Path::new("test_cases/featurecounts/feature_counts_gene_level_stranded"));
}

#[test]
fn test_case_featurecounts_feature_counts_gene_level_unstranded() {
    run_test(std::path::Path::new("test_cases/featurecounts/feature_counts_gene_level_unstranded"));
}

#[test]
fn test_case_htseq_simple_intersection_non_empty() {
    run_test(std::path::Path::new("test_cases/htseq/simple-intersection-non-empty"));
}

#[test]
fn test_case_htseq_simple_intersection_non_empty_non_unique_none_secondary_yes() {
    run_test(std::path::Path::new("test_cases/htseq/simple-intersection-non-empty-non-unique=none-secondary=yes"));
}

#[test]
fn test_case_htseq_simple_intersection_strict() {
    run_test(std::path::Path::new("test_cases/htseq/simple-intersection-strict"));
}

#[test]
fn test_case_htseq_simple_union() {
    run_test(std::path::Path::new("test_cases/htseq/simple-union"));
}

#[test]
fn test_case_order_order_1a() {
    run_test(std::path::Path::new("test_cases/order/order_1a"));
}

#[test]
fn test_case_order_order_1b() {
    run_test(std::path::Path::new("test_cases/order/order_1b"));
}

#[test]
fn test_case_per_region_dedup_per_direction() {
    run_test(std::path::Path::new("test_cases/per_region/dedup_per_direction"));
}

#[test]
fn test_case_per_region_intersecting_regions_aborts() {
    run_test(std::path::Path::new("test_cases/per_region/intersecting_regions_aborts"));
}

#[test]
fn test_case_source_reference_dedup_per_direction() {
    run_test(std::path::Path::new("test_cases/source_reference/dedup_per_direction"));
}

#[test]
fn test_case_source_reference_reference() {
    run_test(std::path::Path::new("test_cases/source_reference/reference"));
}

#[test]
fn test_case_source_reference_reference_forward_only() {
    run_test(std::path::Path::new("test_cases/source_reference/reference_forward_only"));
}

#[test]
fn test_case_source_reference_reference_no_index() {
    run_test(std::path::Path::new("test_cases/source_reference/reference_no_index"));
}

#[test]
fn test_case_source_reference_reference_reverse_only() {
    run_test(std::path::Path::new("test_cases/source_reference/reference_reverse_only"));
}

#[test]
fn test_case_specific_external_umi_thresholder_command() {
    run_test(std::path::Path::new("test_cases/specific/external_umi_thresholder_command"));
}

#[test]
fn test_case_specific_external_umi_thresholder_command_with_barcodes() {
    run_test(std::path::Path::new("test_cases/specific/external_umi_thresholder_command_with_barcodes"));
}

#[test]
fn test_case_specific_filters_min_aligned_len() {
    run_test(std::path::Path::new("test_cases/specific/filters/min_aligned_len"));
}

#[test]
fn test_case_specific_read_clipped_before_start() {
    run_test(std::path::Path::new("test_cases/specific/read_clipped_before_start"));
}

#[test]
fn test_case_starsolo_starsolo_targeted() {
    run_test(std::path::Path::new("test_cases/starsolo/starSolo_targeted"));
}

#[test]
fn test_case_umi_tools_position_based_dedup_cluster() {
    run_test(std::path::Path::new("test_cases/umi-tools/position_based_dedup/cluster"));
}

#[test]
fn test_case_umi_tools_position_based_dedup_directional() {
    run_test(std::path::Path::new("test_cases/umi-tools/position_based_dedup/directional"));
}

#[test]
fn test_case_umi_tools_position_based_dedup_directional_2mm() {
    run_test(std::path::Path::new("test_cases/umi-tools/position_based_dedup/directional_2mm"));
}

#[test]
fn test_case_umi_tools_position_based_dedup_directional_starsolo() {
    run_test(std::path::Path::new("test_cases/umi-tools/position_based_dedup/directional_starsolo"));
}

#[test]
fn test_case_umi_tools_position_based_dedup_percentile() {
    run_test(std::path::Path::new("test_cases/umi-tools/position_based_dedup/percentile"));
}

#[test]
fn test_case_umi_tools_position_based_dedup_unique() {
    run_test(std::path::Path::new("test_cases/umi-tools/position_based_dedup/unique"));
}

#[test]
fn test_case_umi_tools_unique_gene_tag() {
    run_test(std::path::Path::new("test_cases/umi-tools/unique_gene_tag"));
}

#[test]
fn test_case_umi_tools_unique_gene_tag_umi_search_in_name() {
    run_test(std::path::Path::new("test_cases/umi-tools/unique_gene_tag-umi_search_in_name"));
}

#[test]
fn test_case_umi_tools_unique_stranded() {
    run_test(std::path::Path::new("test_cases/umi-tools/unique_stranded"));
}

#[test]
fn test_case_validation_basic_no_gtf_set() {
    run_test(std::path::Path::new("test_cases/validation/basic_no_gtf_set"));
}

#[test]
fn test_case_validation_empty_bam() {
    run_test(std::path::Path::new("test_cases/validation/empty_bam"));
}

#[test]
fn test_case_validation_filter_references_min_length() {
    run_test(std::path::Path::new("test_cases/validation/filter_references_min_length"));
}

#[test]
fn test_case_validation_invalid_bam_tag_0() {
    run_test(std::path::Path::new("test_cases/validation/invalid_bam_tag_0"));
}

#[test]
fn test_case_validation_invalid_bam_tag_1() {
    run_test(std::path::Path::new("test_cases/validation/invalid_bam_tag_1"));
}
