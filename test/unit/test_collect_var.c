#include <stdlib.h>
#include <string.h>

#include "unity.h"

#include "bam_utils.h"
#include "call_var_main.h"
#include "collect_var.h"

cand_var_t *init_cand_vars_based_on_sites(int n_var_sites, var_site_t *var_sites);
int merge_var_sites(const call_var_opt_t *opt, int n_total_var_sites, var_site_t **var_sites, int tid,
                    hts_pos_t reg_beg, hts_pos_t reg_end, int n_digar, digar1_t *digars);
int collect_all_cand_var_sites(const call_var_opt_t *opt, bam_chunk_t *chunk, var_site_t **var_sites);
int merge_var_profile(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_new_vars, cand_var_t *new_vars, int *new_var_cate, read_var_profile_t *new_p);

static cgranges_t *make_read_var_cr(read_var_profile_t *p, int n_reads, int *ordered_read_ids, uint8_t *is_skipped) {
    cgranges_t *read_var_cr = cr_init();

    for (int i = 0; i < n_reads; ++i) {
        int read_i = ordered_read_ids[i];
        if (is_skipped[read_i]) continue;
        if (p[read_i].start_var_idx < 0 || p[read_i].end_var_idx < 0) continue;
        cr_add(read_var_cr, "cr", p[read_i].start_var_idx, p[read_i].end_var_idx + 1, read_i);
    }
    cr_index(read_var_cr);
    return read_var_cr;
}

static void free_merge_chunk(bam_chunk_t *chunk) {
    free_cand_vars(chunk->cand_vars, chunk->n_cand_vars);
    free_read_var_profile(chunk->read_var_profile, chunk->n_reads);
    free(chunk->var_i_to_cate);
    cr_destroy(chunk->read_var_cr);
}

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

static void test_merge_var_profile_interleaves_variants_and_remaps_read_profiles(void) {
    call_var_opt_t opt = {0};
    bam_chunk_t chunk = {0};
    int ordered_read_ids[] = {0, 1, 2};
    uint8_t is_skipped[] = {0, 0, 0};
    var_site_t old_sites[2];
    var_site_t new_sites[2];
    uint8_t old_alt0[] = {'T'};
    uint8_t new_alt0[] = {'A', 'C'};
    uint8_t new_alt1[] = {'G'};
    int *new_var_cate = (int *)malloc(2 * sizeof(int));
    read_var_profile_t *new_p = NULL;

    opt.min_sv_len = 50;
    chunk.n_reads = 3;
    chunk.n_cand_vars = 2;
    chunk.ordered_read_ids = ordered_read_ids;
    chunk.is_skipped = is_skipped;

    old_sites[0] = (var_site_t){0, 100, BAM_CDIFF, 1, 1, old_alt0};
    old_sites[1] = (var_site_t){0, 120, BAM_CDEL, 2, 0, NULL};
    new_sites[0] = (var_site_t){0, 110, BAM_CINS, 0, 2, new_alt0};
    new_sites[1] = (var_site_t){0, 130, BAM_CDIFF, 1, 1, new_alt1};

    chunk.cand_vars = init_cand_vars_based_on_sites(2, old_sites);
    chunk.var_i_to_cate = (int *)malloc(2 * sizeof(int));
    chunk.read_var_profile = init_read_var_profile(3, 2);
    TEST_ASSERT_NOT_NULL(chunk.cand_vars);
    TEST_ASSERT_NOT_NULL(chunk.var_i_to_cate);
    TEST_ASSERT_NOT_NULL(chunk.read_var_profile);
    TEST_ASSERT_NOT_NULL(new_var_cate);

    chunk.var_i_to_cate[0] = 101;
    chunk.var_i_to_cate[1] = 202;
    update_read_var_profile_with_allele(0, 1, 10, chunk.read_var_profile + 0);
    update_read_var_profile_with_allele(1, 0, -1, chunk.read_var_profile + 0);
    update_read_var_profile_with_allele(1, 1, 20, chunk.read_var_profile + 1);
    chunk.read_var_cr = make_read_var_cr(chunk.read_var_profile, chunk.n_reads, chunk.ordered_read_ids, chunk.is_skipped);

    new_p = init_read_var_profile(3, 2);
    TEST_ASSERT_NOT_NULL(new_p);
    new_var_cate[0] = 303;
    new_var_cate[1] = 404;
    update_read_var_profile_with_allele(0, 0, -1, new_p + 0);
    update_read_var_profile_with_allele(0, 1, 22, new_p + 1);
    update_read_var_profile_with_allele(1, 1, 23, new_p + 1);

    TEST_ASSERT_EQUAL_INT(2, merge_var_profile(&opt, &chunk, 2, init_cand_vars_based_on_sites(2, new_sites), new_var_cate, new_p));

    TEST_ASSERT_EQUAL_INT(4, chunk.n_cand_vars);
    TEST_ASSERT_EQUAL_INT(100, chunk.cand_vars[0].pos);
    TEST_ASSERT_EQUAL_INT(110, chunk.cand_vars[1].pos);
    TEST_ASSERT_EQUAL_INT(120, chunk.cand_vars[2].pos);
    TEST_ASSERT_EQUAL_INT(130, chunk.cand_vars[3].pos);
    TEST_ASSERT_EQUAL_INT(101, chunk.var_i_to_cate[0]);
    TEST_ASSERT_EQUAL_INT(303, chunk.var_i_to_cate[1]);
    TEST_ASSERT_EQUAL_INT(202, chunk.var_i_to_cate[2]);
    TEST_ASSERT_EQUAL_INT(404, chunk.var_i_to_cate[3]);

    TEST_ASSERT_EQUAL_INT(0, chunk.read_var_profile[0].start_var_idx);
    TEST_ASSERT_EQUAL_INT(2, chunk.read_var_profile[0].end_var_idx);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[0].alleles[0]);
    TEST_ASSERT_EQUAL_INT(10, chunk.read_var_profile[0].alt_qi[0]);
    TEST_ASSERT_EQUAL_INT(0, chunk.read_var_profile[0].alleles[1]);
    TEST_ASSERT_EQUAL_INT(-1, chunk.read_var_profile[0].alt_qi[1]);
    TEST_ASSERT_EQUAL_INT(0, chunk.read_var_profile[0].alleles[2]);

    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[1].start_var_idx);
    TEST_ASSERT_EQUAL_INT(3, chunk.read_var_profile[1].end_var_idx);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[1].alleles[0]);
    TEST_ASSERT_EQUAL_INT(22, chunk.read_var_profile[1].alt_qi[0]);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[1].alleles[1]);
    TEST_ASSERT_EQUAL_INT(20, chunk.read_var_profile[1].alt_qi[1]);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[1].alleles[2]);
    TEST_ASSERT_EQUAL_INT(23, chunk.read_var_profile[1].alt_qi[2]);

    TEST_ASSERT_EQUAL_INT(-1, chunk.read_var_profile[2].start_var_idx);
    TEST_ASSERT_EQUAL_INT(-2, chunk.read_var_profile[2].end_var_idx);

    TEST_ASSERT_EQUAL_INT64(2, chunk.read_var_cr->n_r);
    TEST_ASSERT_EQUAL_INT(0, cr_start(chunk.read_var_cr, 0));
    TEST_ASSERT_EQUAL_INT(3, cr_end(chunk.read_var_cr, 0));
    TEST_ASSERT_EQUAL_INT(0, cr_label(chunk.read_var_cr, 0));
    TEST_ASSERT_EQUAL_INT(1, cr_start(chunk.read_var_cr, 1));
    TEST_ASSERT_EQUAL_INT(4, cr_end(chunk.read_var_cr, 1));
    TEST_ASSERT_EQUAL_INT(1, cr_label(chunk.read_var_cr, 1));

    free_merge_chunk(&chunk);
}

static void test_merge_var_profile_keeps_old_duplicate_and_skips_skipped_reads(void) {
    call_var_opt_t opt = {0};
    bam_chunk_t chunk = {0};
    int ordered_read_ids[] = {0, 1, 2};
    uint8_t is_skipped[] = {0, 1, 0};
    var_site_t old_sites[1];
    var_site_t new_sites[2];
    uint8_t dup_alt[] = {'T'};
    uint8_t new_alt[] = {'G', 'G'};
    int *new_var_cate = (int *)malloc(2 * sizeof(int));
    read_var_profile_t *new_p = NULL;

    opt.min_sv_len = 50;
    chunk.n_reads = 3;
    chunk.n_cand_vars = 1;
    chunk.ordered_read_ids = ordered_read_ids;
    chunk.is_skipped = is_skipped;

    old_sites[0] = (var_site_t){0, 100, BAM_CDIFF, 1, 1, dup_alt};
    new_sites[0] = (var_site_t){0, 100, BAM_CDIFF, 1, 1, dup_alt};
    new_sites[1] = (var_site_t){0, 105, BAM_CINS, 0, 2, new_alt};

    chunk.cand_vars = init_cand_vars_based_on_sites(1, old_sites);
    chunk.var_i_to_cate = (int *)malloc(sizeof(int));
    chunk.read_var_profile = init_read_var_profile(3, 1);
    TEST_ASSERT_NOT_NULL(chunk.cand_vars);
    TEST_ASSERT_NOT_NULL(chunk.var_i_to_cate);
    TEST_ASSERT_NOT_NULL(chunk.read_var_profile);
    TEST_ASSERT_NOT_NULL(new_var_cate);

    chunk.var_i_to_cate[0] = 11;
    update_read_var_profile_with_allele(0, 1, 10, chunk.read_var_profile + 0);
    update_read_var_profile_with_allele(0, 0, -1, chunk.read_var_profile + 1);
    chunk.read_var_cr = make_read_var_cr(chunk.read_var_profile, chunk.n_reads, chunk.ordered_read_ids, chunk.is_skipped);

    new_p = init_read_var_profile(3, 2);
    TEST_ASSERT_NOT_NULL(new_p);
    new_var_cate[0] = 22;
    new_var_cate[1] = 33;
    update_read_var_profile_with_allele(0, 0, -1, new_p + 0);
    update_read_var_profile_with_allele(1, 1, 15, new_p + 0);
    update_read_var_profile_with_allele(0, 1, 17, new_p + 1);
    update_read_var_profile_with_allele(1, 1, 18, new_p + 1);
    update_read_var_profile_with_allele(1, 0, -1, new_p + 2);

    TEST_ASSERT_EQUAL_INT(2, merge_var_profile(&opt, &chunk, 2, init_cand_vars_based_on_sites(2, new_sites), new_var_cate, new_p));

    TEST_ASSERT_EQUAL_INT(2, chunk.n_cand_vars);
    TEST_ASSERT_EQUAL_INT(100, chunk.cand_vars[0].pos);
    TEST_ASSERT_EQUAL_INT(105, chunk.cand_vars[1].pos);
    TEST_ASSERT_EQUAL_INT(11, chunk.var_i_to_cate[0]);
    TEST_ASSERT_EQUAL_INT(33, chunk.var_i_to_cate[1]);

    TEST_ASSERT_EQUAL_INT(0, chunk.read_var_profile[0].start_var_idx);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[0].end_var_idx);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[0].alleles[0]);
    TEST_ASSERT_EQUAL_INT(10, chunk.read_var_profile[0].alt_qi[0]);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[0].alleles[1]);
    TEST_ASSERT_EQUAL_INT(15, chunk.read_var_profile[0].alt_qi[1]);

    TEST_ASSERT_EQUAL_INT(-1, chunk.read_var_profile[1].start_var_idx);
    TEST_ASSERT_EQUAL_INT(-2, chunk.read_var_profile[1].end_var_idx);

    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[2].start_var_idx);
    TEST_ASSERT_EQUAL_INT(1, chunk.read_var_profile[2].end_var_idx);
    TEST_ASSERT_EQUAL_INT(0, chunk.read_var_profile[2].alleles[0]);
    TEST_ASSERT_EQUAL_INT(-1, chunk.read_var_profile[2].alt_qi[0]);

    TEST_ASSERT_EQUAL_INT64(2, chunk.read_var_cr->n_r);
    TEST_ASSERT_EQUAL_INT(0, cr_label(chunk.read_var_cr, 0));
    TEST_ASSERT_EQUAL_INT(2, cr_label(chunk.read_var_cr, 1));

    free_merge_chunk(&chunk);
}

void collect_var_suite(void) {
    RUN_TEST(test_merge_var_sites_filters_and_keeps_sorted_unique_sites);
    RUN_TEST(test_merge_var_sites_fuzzy_merges_large_insertions);
    RUN_TEST(test_collect_all_cand_var_sites_batches_across_reads);
    RUN_TEST(test_merge_var_profile_interleaves_variants_and_remaps_read_profiles);
    RUN_TEST(test_merge_var_profile_keeps_old_duplicate_and_skips_skipped_reads);
}
