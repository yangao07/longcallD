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
int exact_comp_cand_var(const call_var_opt_t *opt, cand_var_t *var1, cand_var_t *var2);
float var_noisy_reads_ratio(bam_chunk_t *chunk, hts_pos_t var_start, hts_pos_t var_end);
void free_cand_vars1(cand_var_t *cand_vars);

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

static read_var_profile_t *dup_read_var_profile(const read_var_profile_t *src, int n_reads, int n_total_vars) {
    read_var_profile_t *dst = init_read_var_profile(n_reads, n_total_vars);

    TEST_ASSERT_NOT_NULL(dst);
    for (int i = 0; i < n_reads; ++i) {
        dst[i].read_id = src[i].read_id;
        dst[i].start_var_idx = src[i].start_var_idx;
        dst[i].end_var_idx = src[i].end_var_idx;
        memcpy(dst[i].alleles, src[i].alleles, (size_t)n_total_vars * sizeof(int));
        memcpy(dst[i].alt_qi, src[i].alt_qi, (size_t)n_total_vars * sizeof(int));
    }
    return dst;
}

static int *dup_int_array(const int *src, int n) {
    int *dst = (int *)malloc((size_t)n * sizeof(int));

    TEST_ASSERT_NOT_NULL(dst);
    memcpy(dst, src, (size_t)n * sizeof(int));
    return dst;
}

static uint32_t next_lcg(uint32_t *state) {
    *state = (*state * 1664525u) + 1013904223u;
    return *state;
}

static int choose_allele(uint32_t *state) {
    static const int alleles[] = {-1, 0, 1};
    return alleles[next_lcg(state) % 3];
}

static int merge_var_profile_reference(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_new_vars, cand_var_t *new_vars, int *new_var_cate, read_var_profile_t *new_p) {
    if (n_new_vars <= 0) return 0;
    cand_var_t *old_vars = chunk->cand_vars;
    read_var_profile_t *old_p = chunk->read_var_profile;
    cand_var_t *merged_vars = (cand_var_t *)malloc((chunk->n_cand_vars + n_new_vars) * sizeof(cand_var_t));
    read_var_profile_t *merged_p = init_read_var_profile(chunk->n_reads, chunk->n_cand_vars + n_new_vars);
    int *merged_var_i_to_cate = (int*)malloc((chunk->n_cand_vars + n_new_vars) * sizeof(int));
    int old_var_i = 0, new_var_i = 0, merged_var_i = 0;

    TEST_ASSERT_NOT_NULL(merged_vars);
    TEST_ASSERT_NOT_NULL(merged_p);
    TEST_ASSERT_NOT_NULL(merged_var_i_to_cate);

    for (; old_var_i < chunk->n_cand_vars && new_var_i < n_new_vars; ) {
        int ret = exact_comp_cand_var(opt, old_vars+old_var_i, new_vars+new_var_i);
        if (ret < 0) {
            for (int i = 0; i < chunk->n_reads; ++i) {
                int read_i = chunk->ordered_read_ids[i];
                if (chunk->is_skipped[read_i]) continue;
                read_var_profile_t *old_p1 = old_p + read_i;
                read_var_profile_t *merged_p1 = merged_p + read_i;
                if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], old_p1->alt_qi[old_var_i-old_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
            merged_vars[merged_var_i++] = old_vars[old_var_i++];
        } else if (ret > 0) {
            for (int i = 0; i < chunk->n_reads; ++i) {
                int read_i = chunk->ordered_read_ids[i];
                if (chunk->is_skipped[read_i]) continue;
                read_var_profile_t *new_p1 = new_p + read_i;
                read_var_profile_t *merged_p1 = merged_p + read_i;
                if (new_p1->start_var_idx > new_var_i || new_p1->end_var_idx < new_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, new_p1->alleles[new_var_i-new_p1->start_var_idx], new_p1->alt_qi[new_var_i-new_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = new_var_cate[new_var_i];
            merged_vars[merged_var_i++] = new_vars[new_var_i++];
        } else {
            for (int i = 0; i < chunk->n_reads; ++i) {
                int read_i = chunk->ordered_read_ids[i];
                if (chunk->is_skipped[read_i]) continue;
                read_var_profile_t *old_p1 = old_p + read_i;
                read_var_profile_t *merged_p1 = merged_p + read_i;
                if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], old_p1->alt_qi[old_var_i-old_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
            merged_vars[merged_var_i++] = old_vars[old_var_i++];
            free_cand_vars1(new_vars+new_var_i);
            new_var_i++;
        }
    }
    for (; old_var_i < chunk->n_cand_vars; ++old_var_i) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            int read_i = chunk->ordered_read_ids[i];
            if (chunk->is_skipped[read_i]) continue;
            read_var_profile_t *old_p1 = old_p + read_i;
            read_var_profile_t *merged_p1 = merged_p + read_i;
            if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
            update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], old_p1->alt_qi[old_var_i-old_p1->start_var_idx], merged_p1);
        }
        merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
        merged_vars[merged_var_i++] = old_vars[old_var_i];
    }
    for (; new_var_i < n_new_vars; ++new_var_i) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            int read_i = chunk->ordered_read_ids[i];
            if (chunk->is_skipped[read_i]) continue;
            read_var_profile_t *new_p1 = new_p + read_i;
            read_var_profile_t *merged_p1 = merged_p + read_i;
            if (new_p1->start_var_idx > new_var_i || new_p1->end_var_idx < new_var_i) continue;
            update_read_var_profile_with_allele(merged_var_i, new_p1->alleles[new_var_i-new_p1->start_var_idx], new_p1->alt_qi[new_var_i-new_p1->start_var_idx], merged_p1);
        }
        merged_var_i_to_cate[merged_var_i] = new_var_cate[new_var_i];
        merged_vars[merged_var_i++] = new_vars[new_var_i];
    }

    cgranges_t *merged_read_var_cr = cr_init();
    TEST_ASSERT_NOT_NULL(merged_read_var_cr);
    for (int i = 0; i < chunk->n_reads; ++i) {
        int read_i = chunk->ordered_read_ids[i];
        if (chunk->is_skipped[read_i]) continue;
        if (merged_p[read_i].start_var_idx < 0 || merged_p[read_i].end_var_idx < 0) continue;
        cr_add(merged_read_var_cr, "cr", merged_p[read_i].start_var_idx, merged_p[read_i].end_var_idx+1, read_i);
    }
    cr_index(merged_read_var_cr);

    if (chunk->n_cand_vars > 0) free_read_var_profile(chunk->read_var_profile, chunk->n_reads);
    free_read_var_profile(new_p, chunk->n_reads);
    free(chunk->var_i_to_cate);
    free(chunk->cand_vars);
    free(new_vars);
    free(new_var_cate);
    cr_destroy(chunk->read_var_cr);

    chunk->read_var_profile = merged_p;
    chunk->cand_vars = merged_vars;
    chunk->n_cand_vars = merged_var_i;
    chunk->var_i_to_cate = merged_var_i_to_cate;
    chunk->read_var_cr = merged_read_var_cr;
    return new_var_i;
}

static void assert_merge_chunks_equal(const bam_chunk_t *expected, const bam_chunk_t *actual) {
    TEST_ASSERT_EQUAL_INT(expected->n_cand_vars, actual->n_cand_vars);
    for (int i = 0; i < expected->n_cand_vars; ++i) {
        TEST_ASSERT_EQUAL_INT64(expected->cand_vars[i].pos, actual->cand_vars[i].pos);
        TEST_ASSERT_EQUAL_INT(expected->cand_vars[i].var_type, actual->cand_vars[i].var_type);
        TEST_ASSERT_EQUAL_INT(expected->cand_vars[i].ref_len, actual->cand_vars[i].ref_len);
        TEST_ASSERT_EQUAL_INT(expected->cand_vars[i].alt_len, actual->cand_vars[i].alt_len);
        TEST_ASSERT_EQUAL_INT(expected->var_i_to_cate[i], actual->var_i_to_cate[i]);
        if (expected->cand_vars[i].alt_len > 0) {
            TEST_ASSERT_EQUAL_MEMORY(expected->cand_vars[i].alt_seq, actual->cand_vars[i].alt_seq, (size_t)expected->cand_vars[i].alt_len);
        } else {
            TEST_ASSERT_NULL(actual->cand_vars[i].alt_seq);
        }
    }

    for (int read_i = 0; read_i < expected->n_reads; ++read_i) {
        TEST_ASSERT_EQUAL_INT(expected->read_var_profile[read_i].start_var_idx, actual->read_var_profile[read_i].start_var_idx);
        TEST_ASSERT_EQUAL_INT(expected->read_var_profile[read_i].end_var_idx, actual->read_var_profile[read_i].end_var_idx);
        TEST_ASSERT_EQUAL_MEMORY(expected->read_var_profile[read_i].alleles, actual->read_var_profile[read_i].alleles,
                                 (size_t)expected->n_cand_vars * sizeof(int));
        TEST_ASSERT_EQUAL_MEMORY(expected->read_var_profile[read_i].alt_qi, actual->read_var_profile[read_i].alt_qi,
                                 (size_t)expected->n_cand_vars * sizeof(int));
    }

    TEST_ASSERT_EQUAL_INT64(expected->read_var_cr->n_r, actual->read_var_cr->n_r);
    for (int i = 0; i < expected->read_var_cr->n_r; ++i) {
        TEST_ASSERT_EQUAL_INT(cr_start(expected->read_var_cr, i), cr_start(actual->read_var_cr, i));
        TEST_ASSERT_EQUAL_INT(cr_end(expected->read_var_cr, i), cr_end(actual->read_var_cr, i));
        TEST_ASSERT_EQUAL_INT(cr_label(expected->read_var_cr, i), cr_label(actual->read_var_cr, i));
    }
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

static void test_merge_var_profile_matches_reference_on_sparse_profiles(void) {
    call_var_opt_t opt = {0};
    uint32_t rng = 7;
    int ordered_read_ids[] = {0, 1, 2, 3, 4, 5};
    uint8_t is_skipped[6];

    opt.min_sv_len = 50;

    for (int iter = 0; iter < 32; ++iter) {
        enum { MAX_SITES = 8, N_READS = 6 };
        bam_chunk_t expected = {0}, actual = {0};
        var_site_t old_sites[MAX_SITES];
        var_site_t new_sites[MAX_SITES];
        int old_var_cate[MAX_SITES];
        int new_var_cate[MAX_SITES];
        uint8_t old_alt_buf[MAX_SITES][2];
        uint8_t new_alt_buf[MAX_SITES][2];
        int n_old_vars = 0, n_new_vars = 0;
        read_var_profile_t *old_profiles = NULL, *new_profiles = NULL;

        for (int read_i = 0; read_i < N_READS; ++read_i) {
            is_skipped[read_i] = (uint8_t)(next_lcg(&rng) % 5 == 0);
        }

        for (int site_i = 0; site_i < MAX_SITES; ++site_i) {
            int include_old = (next_lcg(&rng) & 1u) != 0;
            int include_new = (next_lcg(&rng) & 1u) != 0;
            int pos = 100 + site_i * 10;

            if (!include_old && !include_new) include_old = 1;
            if (site_i == 0 && !include_new) include_new = 1;

            if (include_old) {
                old_alt_buf[n_old_vars][0] = (uint8_t)("TGCA"[site_i % 4] - 'A');
                old_alt_buf[n_old_vars][1] = (uint8_t)("ACGT"[(site_i + 1) % 4] - 'A');
                old_sites[n_old_vars].tid = 0;
                old_sites[n_old_vars].pos = pos;
                old_sites[n_old_vars].var_type = (site_i % 3 == 0 ? BAM_CDIFF : (site_i % 3 == 1 ? BAM_CINS : BAM_CDEL));
                old_sites[n_old_vars].ref_len = (old_sites[n_old_vars].var_type == BAM_CDEL ? 2 : (old_sites[n_old_vars].var_type == BAM_CDIFF ? 1 : 0));
                old_sites[n_old_vars].alt_len = (old_sites[n_old_vars].var_type == BAM_CINS ? 2 : (old_sites[n_old_vars].var_type == BAM_CDIFF ? 1 : 0));
                old_sites[n_old_vars].alt_seq = (old_sites[n_old_vars].alt_len > 0 ? old_alt_buf[n_old_vars] : NULL);
                old_var_cate[n_old_vars] = 1000 + site_i;
                n_old_vars++;
            }

            if (include_new) {
                if (include_old && (next_lcg(&rng) % 3 == 0)) {
                    new_sites[n_new_vars] = old_sites[n_old_vars-1];
                    if (new_sites[n_new_vars].alt_len > 0) {
                        memcpy(new_alt_buf[n_new_vars], old_alt_buf[n_old_vars-1], (size_t)new_sites[n_new_vars].alt_len);
                        new_sites[n_new_vars].alt_seq = new_alt_buf[n_new_vars];
                    }
                } else {
                    new_alt_buf[n_new_vars][0] = (uint8_t)("CATG"[(site_i + 2) % 4] - 'A');
                    new_alt_buf[n_new_vars][1] = (uint8_t)("GTAC"[(site_i + 3) % 4] - 'A');
                    new_sites[n_new_vars].tid = 0;
                    new_sites[n_new_vars].pos = pos;
                    new_sites[n_new_vars].var_type = (site_i % 2 == 0 ? BAM_CDIFF : BAM_CINS);
                    new_sites[n_new_vars].ref_len = (new_sites[n_new_vars].var_type == BAM_CDIFF ? 1 : 0);
                    new_sites[n_new_vars].alt_len = (new_sites[n_new_vars].var_type == BAM_CDIFF ? 1 : 2);
                    new_sites[n_new_vars].alt_seq = new_alt_buf[n_new_vars];
                }
                new_var_cate[n_new_vars] = 2000 + site_i;
                n_new_vars++;
            }
        }

        old_profiles = init_read_var_profile(N_READS, n_old_vars);
        new_profiles = init_read_var_profile(N_READS, n_new_vars);
        TEST_ASSERT_NOT_NULL(old_profiles);
        TEST_ASSERT_NOT_NULL(new_profiles);

        for (int read_i = 0; read_i < N_READS; ++read_i) {
            if (!is_skipped[read_i]) {
                for (int var_i = 0; var_i < n_old_vars; ++var_i) {
                    if ((next_lcg(&rng) % 3) != 0) {
                        update_read_var_profile_with_allele(var_i, choose_allele(&rng), (int)(next_lcg(&rng) % 50), old_profiles + read_i);
                    }
                }
                for (int var_i = 0; var_i < n_new_vars; ++var_i) {
                    if ((next_lcg(&rng) % 3) != 0) {
                        update_read_var_profile_with_allele(var_i, choose_allele(&rng), (int)(next_lcg(&rng) % 50), new_profiles + read_i);
                    }
                }
            }
        }

        expected.n_reads = N_READS;
        expected.n_cand_vars = n_old_vars;
        expected.ordered_read_ids = ordered_read_ids;
        expected.is_skipped = is_skipped;
        expected.cand_vars = init_cand_vars_based_on_sites(n_old_vars, old_sites);
        expected.var_i_to_cate = dup_int_array(old_var_cate, n_old_vars);
        expected.read_var_profile = dup_read_var_profile(old_profiles, N_READS, n_old_vars);
        expected.read_var_cr = make_read_var_cr(expected.read_var_profile, N_READS, expected.ordered_read_ids, expected.is_skipped);

        actual.n_reads = N_READS;
        actual.n_cand_vars = n_old_vars;
        actual.ordered_read_ids = ordered_read_ids;
        actual.is_skipped = is_skipped;
        actual.cand_vars = init_cand_vars_based_on_sites(n_old_vars, old_sites);
        actual.var_i_to_cate = dup_int_array(old_var_cate, n_old_vars);
        actual.read_var_profile = dup_read_var_profile(old_profiles, N_READS, n_old_vars);
        actual.read_var_cr = make_read_var_cr(actual.read_var_profile, N_READS, actual.ordered_read_ids, actual.is_skipped);

        TEST_ASSERT_EQUAL_INT(n_new_vars, merge_var_profile_reference(&opt, &expected, n_new_vars,
                               init_cand_vars_based_on_sites(n_new_vars, new_sites), dup_int_array(new_var_cate, n_new_vars),
                               dup_read_var_profile(new_profiles, N_READS, n_new_vars)));
        TEST_ASSERT_EQUAL_INT(n_new_vars, merge_var_profile(&opt, &actual, n_new_vars,
                               init_cand_vars_based_on_sites(n_new_vars, new_sites), dup_int_array(new_var_cate, n_new_vars),
                               dup_read_var_profile(new_profiles, N_READS, n_new_vars)));

        assert_merge_chunks_equal(&expected, &actual);

        free_merge_chunk(&expected);
        free_merge_chunk(&actual);
        free_read_var_profile(old_profiles, N_READS);
        free_read_var_profile(new_profiles, N_READS);
    }
}

static void test_var_noisy_reads_ratio_uses_interval_cache_and_dedupes_reads(void) {
    bam_chunk_t chunk;
    int ordered_read_ids[] = {0, 1, 2};
    uint8_t is_skipped[] = {0, 0, 1};
    digar1_t *read0_digars = (digar1_t *)calloc(2, sizeof(digar1_t));
    digar1_t *read1_digars = (digar1_t *)calloc(1, sizeof(digar1_t));
    digar1_t *read2_digars = (digar1_t *)calloc(1, sizeof(digar1_t));

    memset(&chunk, 0, sizeof(chunk));
    TEST_ASSERT_NOT_NULL(read0_digars);
    TEST_ASSERT_NOT_NULL(read1_digars);
    TEST_ASSERT_NOT_NULL(read2_digars);

    chunk.n_reads = 3;
    chunk.m_reads = 3;
    chunk.ordered_read_ids = ordered_read_ids;
    chunk.is_skipped = is_skipped;
    chunk.digars = (digar_t *)calloc(3, sizeof(digar_t));
    chunk.tname = (char *)"chr20";
    TEST_ASSERT_NOT_NULL(chunk.digars);

    chunk.digars[0].beg = 100;
    chunk.digars[0].end = 120;
    chunk.digars[0].n_digar = chunk.digars[0].m_digar = 2;
    chunk.digars[0].digars = read0_digars;
    read0_digars[0] = (digar1_t){105, BAM_CDIFF, 1, 0, NULL, 0};
    read0_digars[1] = (digar1_t){110, BAM_CDIFF, 1, 5, NULL, 0};

    chunk.digars[1].beg = 100;
    chunk.digars[1].end = 120;
    chunk.digars[1].n_digar = chunk.digars[1].m_digar = 1;
    chunk.digars[1].digars = read1_digars;
    read1_digars[0] = (digar1_t){115, BAM_CDIFF, 1, 0, NULL, 0};

    chunk.digars[2].beg = 130;
    chunk.digars[2].end = 150;
    chunk.digars[2].n_digar = chunk.digars[2].m_digar = 1;
    chunk.digars[2].digars = read2_digars;
    read2_digars[0] = (digar1_t){135, BAM_CDIFF, 1, 0, NULL, 0};

    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.5f, var_noisy_reads_ratio(&chunk, 105, 105));
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.5f, var_noisy_reads_ratio(&chunk, 105, 110));
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 1.0f, var_noisy_reads_ratio(&chunk, 105, 116));
    TEST_ASSERT_FLOAT_WITHIN(1e-6f, 0.0f, var_noisy_reads_ratio(&chunk, 135, 135));

    cr_destroy(chunk.var_noisy_read_cov_cr);
    cr_destroy(chunk.var_noisy_read_err_cr);
    free(chunk.var_noisy_read_marks);
    free_digar1(read0_digars, 2);
    free_digar1(read1_digars, 1);
    free_digar1(read2_digars, 1);
    free(chunk.digars);
}

void collect_var_suite(void) {
    RUN_TEST(test_merge_var_sites_filters_and_keeps_sorted_unique_sites);
    RUN_TEST(test_merge_var_sites_fuzzy_merges_large_insertions);
    RUN_TEST(test_collect_all_cand_var_sites_batches_across_reads);
    RUN_TEST(test_merge_var_profile_interleaves_variants_and_remaps_read_profiles);
    RUN_TEST(test_merge_var_profile_keeps_old_duplicate_and_skips_skipped_reads);
    RUN_TEST(test_merge_var_profile_matches_reference_on_sparse_profiles);
    RUN_TEST(test_var_noisy_reads_ratio_uses_interval_cache_and_dedupes_reads);
}
