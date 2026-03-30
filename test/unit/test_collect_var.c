#include <stdlib.h>
#include <string.h>

#include "unity.h"

#include "bam_utils.h"
#include "call_var_main.h"
#include "collect_var.h"

int merge_var_sites(const call_var_opt_t *opt, int n_total_var_sites, var_site_t **var_sites, int tid,
                    hts_pos_t reg_beg, hts_pos_t reg_end, int n_digar, digar1_t *digars);
int collect_all_cand_var_sites(const call_var_opt_t *opt, bam_chunk_t *chunk, var_site_t **var_sites);

static void test_merge_var_sites_filters_and_keeps_sorted_unique_sites(void) {
    call_var_opt_t opt = {0};
    opt.min_sv_len = 50;

    uint8_t snp_alt[] = {'T'};
    var_site_t *var_sites = (var_site_t *)malloc(2 * sizeof(var_site_t));
    digar1_t digars[5];
    uint8_t below_reg_alt[] = {'G'};
    uint8_t duplicate_alt[] = {'T'};
    uint8_t ins_alt[] = {'A', 'C'};
    uint8_t low_qual_alt[] = {'C'};

    TEST_ASSERT_NOT_NULL(var_sites);
    var_sites[0] = (var_site_t){7, 100, BAM_CDIFF, 1, 1, snp_alt};
    var_sites[1] = (var_site_t){7, 110, BAM_CDEL, 2, 0, NULL};

    digars[0] = (digar1_t){99, BAM_CDIFF, 1, 0, below_reg_alt, 0};
    digars[1] = (digar1_t){100, BAM_CDIFF, 1, 1, duplicate_alt, 0};
    digars[2] = (digar1_t){103, BAM_CINS, 2, 2, ins_alt, 0};
    digars[3] = (digar1_t){108, BAM_CDIFF, 1, 4, low_qual_alt, 1};
    digars[4] = (digar1_t){120, BAM_CDEL, 1, 5, NULL, 0};

    int n_sites = merge_var_sites(&opt, 2, &var_sites, 7, 100, 115, 5, digars);

    TEST_ASSERT_EQUAL_INT(3, n_sites);
    TEST_ASSERT_EQUAL_INT(100, var_sites[0].pos);
    TEST_ASSERT_EQUAL_INT(BAM_CDIFF, var_sites[0].var_type);
    TEST_ASSERT_EQUAL_INT(1, var_sites[0].alt_len);
    TEST_ASSERT_EQUAL_MEMORY(snp_alt, var_sites[0].alt_seq, 1);

    TEST_ASSERT_EQUAL_INT(103, var_sites[1].pos);
    TEST_ASSERT_EQUAL_INT(BAM_CINS, var_sites[1].var_type);
    TEST_ASSERT_EQUAL_INT(0, var_sites[1].ref_len);
    TEST_ASSERT_EQUAL_INT(2, var_sites[1].alt_len);
    TEST_ASSERT_EQUAL_MEMORY(ins_alt, var_sites[1].alt_seq, 2);

    TEST_ASSERT_EQUAL_INT(110, var_sites[2].pos);
    TEST_ASSERT_EQUAL_INT(BAM_CDEL, var_sites[2].var_type);
    TEST_ASSERT_EQUAL_INT(2, var_sites[2].ref_len);
    TEST_ASSERT_EQUAL_INT(0, var_sites[2].alt_len);
    TEST_ASSERT_NULL(var_sites[2].alt_seq);

    free(var_sites);
}

static void test_merge_var_sites_fuzzy_merges_large_insertions(void) {
    call_var_opt_t opt = {0};
    opt.min_sv_len = 50;

    uint8_t existing_alt[100];
    uint8_t fuzzy_alt[90];
    var_site_t *var_sites = (var_site_t *)malloc(sizeof(var_site_t));
    digar1_t digars[1];

    memset(existing_alt, 'A', sizeof(existing_alt));
    memset(fuzzy_alt, 'C', sizeof(fuzzy_alt));
    TEST_ASSERT_NOT_NULL(var_sites);

    var_sites[0] = (var_site_t){5, 200, BAM_CINS, 0, (int)sizeof(existing_alt), existing_alt};
    digars[0] = (digar1_t){200, BAM_CINS, (int)sizeof(fuzzy_alt), 0, fuzzy_alt, 0};

    int n_sites = merge_var_sites(&opt, 1, &var_sites, 5, -1, -1, 1, digars);

    TEST_ASSERT_EQUAL_INT(1, n_sites);
    TEST_ASSERT_EQUAL_INT(200, var_sites[0].pos);
    TEST_ASSERT_EQUAL_INT(BAM_CINS, var_sites[0].var_type);
    TEST_ASSERT_EQUAL_INT((int)sizeof(existing_alt), var_sites[0].alt_len);
    TEST_ASSERT_EQUAL_PTR(existing_alt, var_sites[0].alt_seq);

    free(var_sites);
}

static void test_collect_all_cand_var_sites_batches_across_reads(void) {
    call_var_opt_t opt = {0};
    bam_chunk_t chunk = {0};
    int ordered_read_ids[] = {1, 0, 2};
    uint8_t is_skipped[] = {0, 0, 1};
    digar_t digars[3] = {0};
    digar1_t read0_digars[3];
    digar1_t read1_digars[3];
    digar1_t read2_digars[1];
    uint8_t snp_alt[] = {'T'};
    uint8_t ins_alt[] = {'A', 'C'};
    uint8_t skipped_alt[] = {'G'};
    var_site_t *var_sites = NULL;

    opt.min_sv_len = 50;
    chunk.n_reads = 3;
    chunk.tid = 20;
    chunk.reg_beg = 100;
    chunk.reg_end = 115;
    chunk.ordered_read_ids = ordered_read_ids;
    chunk.is_skipped = is_skipped;
    chunk.digars = digars;

    read0_digars[0] = (digar1_t){100, BAM_CDIFF, 1, 0, snp_alt, 0};
    read0_digars[1] = (digar1_t){103, BAM_CINS, 2, 1, ins_alt, 0};
    read0_digars[2] = (digar1_t){118, BAM_CDEL, 1, 3, NULL, 0};
    read1_digars[0] = (digar1_t){99, BAM_CDIFF, 1, 0, skipped_alt, 0};
    read1_digars[1] = (digar1_t){100, BAM_CDIFF, 1, 1, snp_alt, 0};
    read1_digars[2] = (digar1_t){106, BAM_CDIFF, 1, 2, skipped_alt, 1};
    read2_digars[0] = (digar1_t){111, BAM_CDEL, 2, 0, NULL, 0};

    digars[0].n_digar = 3;
    digars[0].digars = read0_digars;
    digars[1].n_digar = 3;
    digars[1].digars = read1_digars;
    digars[2].n_digar = 1;
    digars[2].digars = read2_digars;

    int n_sites = collect_all_cand_var_sites(&opt, &chunk, &var_sites);

    TEST_ASSERT_EQUAL_INT(2, n_sites);
    TEST_ASSERT_NOT_NULL(var_sites);
    TEST_ASSERT_EQUAL_INT(100, var_sites[0].pos);
    TEST_ASSERT_EQUAL_INT(BAM_CDIFF, var_sites[0].var_type);
    TEST_ASSERT_EQUAL_INT(103, var_sites[1].pos);
    TEST_ASSERT_EQUAL_INT(BAM_CINS, var_sites[1].var_type);

    free(var_sites);
}

void collect_var_suite(void) {
    RUN_TEST(test_merge_var_sites_filters_and_keeps_sorted_unique_sites);
    RUN_TEST(test_merge_var_sites_fuzzy_merges_large_insertions);
    RUN_TEST(test_collect_all_cand_var_sites_batches_across_reads);
}
