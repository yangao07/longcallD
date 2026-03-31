#include <stdlib.h>
#include <string.h>

#include "unity.h"

#include "bam_utils.h"
#include "call_var_main.h"

int get_var_start(cand_var_t *var_sites, int cur_site_i, int n_total_pos, hts_pos_t start);
int get_var_site_start(var_site_t *var_sites, int cur_site_i, int n_total_pos, hts_pos_t start);

typedef int (*collect_digar_fn)(bam_chunk_t *chunk, int read_i, const call_var_opt_t *opt, digar_t *digar);

static bam1_t *make_test_insertion_read(int with_md_tag) {
    static const uint32_t cigar[] = {
        bam_cigar_gen(1, BAM_CMATCH),
        bam_cigar_gen(2, BAM_CINS),
        bam_cigar_gen(1, BAM_CMATCH),
    };
    static const char seq[] = "ATGC";
    static const char qual[] = {30, 5, 30, 30};
    bam1_t *read = bam_init1();

    TEST_ASSERT_NOT_NULL(read);
    TEST_ASSERT_GREATER_OR_EQUAL_INT(0, bam_set1(read, sizeof("read1"), "read1", 0, 0, 0, 60,
                                                 sizeof(cigar) / sizeof(cigar[0]), cigar,
                                                 -1, -1, 0, 4, seq, qual, 0));
    if (with_md_tag) {
        static const char md[] = "2";
        TEST_ASSERT_EQUAL_INT(0, bam_aux_append(read, "MD", 'Z', sizeof(md), (uint8_t *)md));
    }
    return read;
}

static void init_test_chunk(bam_chunk_t *chunk, bam1_t *read, const char *ref_seq) {
    memset(chunk, 0, sizeof(*chunk));
    chunk->reads = (bam1_t **)malloc(sizeof(bam1_t *));
    chunk->qual_counts = (int *)calloc(256, sizeof(int));
    chunk->is_ont_palindrome = (uint8_t *)calloc(1, sizeof(uint8_t));
    chunk->chunk_noisy_regs = cr_init();

    TEST_ASSERT_NOT_NULL(chunk->reads);
    TEST_ASSERT_NOT_NULL(chunk->qual_counts);
    TEST_ASSERT_NOT_NULL(chunk->is_ont_palindrome);
    TEST_ASSERT_NOT_NULL(chunk->chunk_noisy_regs);

    chunk->reads[0] = read;
    chunk->n_reads = 1;
    chunk->m_reads = 1;
    chunk->reg_beg = 1;
    chunk->reg_end = 2;
    chunk->ref_beg = 1;
    chunk->ref_end = 2;
    chunk->whole_ref_len = 2;
    chunk->tname = (char *)"chrTest";
    if (ref_seq != NULL) {
        size_t ref_seq_len = strlen(ref_seq) + 1;
        chunk->ref_seq = (char *)malloc(ref_seq_len);
        TEST_ASSERT_NOT_NULL(chunk->ref_seq);
        memcpy(chunk->ref_seq, ref_seq, ref_seq_len);
    }
}

static void free_test_chunk(bam_chunk_t *chunk) {
    if (chunk->reads != NULL) {
        if (chunk->reads[0] != NULL) bam_destroy1(chunk->reads[0]);
        free(chunk->reads);
    }
    free(chunk->qual_counts);
    free(chunk->is_ont_palindrome);
    free(chunk->ref_seq);
    if (chunk->chunk_noisy_regs != NULL) cr_destroy(chunk->chunk_noisy_regs);
}

static void free_test_digar(digar_t *digar) {
    free_digar1(digar->digars, digar->n_digar);
    free(digar->bseq);
    free(digar->qual);
    if (digar->noisy_regs != NULL) cr_destroy(digar->noisy_regs);
}

static void assert_insertion_quality_uses_all_inserted_bases(collect_digar_fn collect_fn, int with_md_tag) {
    call_var_opt_t opt = {0};
    bam_chunk_t chunk;
    digar_t digar = {0};
    bam1_t *read = make_test_insertion_read(with_md_tag);
    digar1_t *ins = NULL;

    opt.min_bq = 20;
    opt.noisy_reg_max_xgaps = 16;
    opt.noisy_reg_slide_win = 10;
    opt.max_noisy_frac_per_read = 2.0;
    opt.max_var_ratio_per_read = 2.0;
    opt.end_clip_reg = 1000;
    init_test_chunk(&chunk, read, "AC");

    TEST_ASSERT_EQUAL_INT(0, collect_fn(&chunk, 0, &opt, &digar));
    TEST_ASSERT_EQUAL_INT(3, digar.n_digar);
    TEST_ASSERT_EQUAL_INT(BAM_CEQUAL, digar.digars[0].type);
    TEST_ASSERT_EQUAL_INT(BAM_CINS, digar.digars[1].type);
    TEST_ASSERT_EQUAL_INT(BAM_CEQUAL, digar.digars[2].type);

    ins = &digar.digars[1];
    TEST_ASSERT_EQUAL_INT(2, ins->len);
    TEST_ASSERT_EQUAL_INT(1, ins->qi);
    TEST_ASSERT_EQUAL_INT(0, ins->is_low_qual);

    free_test_digar(&digar);
    free_test_chunk(&chunk);
}

static void test_init_read_var_profile_initializes_each_read(void) {
    const int n_reads = 3;
    const int n_total_vars = 4;
    read_var_profile_t *p = init_read_var_profile(n_reads, n_total_vars);

    TEST_ASSERT_NOT_NULL(p);
    for (int i = 0; i < n_reads; ++i) {
        TEST_ASSERT_EQUAL_INT(i, p[i].read_id);
        TEST_ASSERT_EQUAL_INT(-1, p[i].start_var_idx);
        TEST_ASSERT_EQUAL_INT(-2, p[i].end_var_idx);
        TEST_ASSERT_NOT_NULL(p[i].alleles);
        TEST_ASSERT_NOT_NULL(p[i].alt_qi);
        for (int j = 0; j < n_total_vars; ++j) {
            TEST_ASSERT_EQUAL_INT(-1, p[i].alleles[j]);
            TEST_ASSERT_EQUAL_INT(-1, p[i].alt_qi[j]);
        }
    }

    free_read_var_profile(p, n_reads);
}

static void test_init_read_var_profile_with_ids_preserves_input_ids(void) {
    int read_ids[] = {7, 3};
    read_var_profile_t *p = init_read_var_profile_with_ids(2, read_ids, 3);

    TEST_ASSERT_NOT_NULL(p);
    TEST_ASSERT_EQUAL_INT(7, p[0].read_id);
    TEST_ASSERT_EQUAL_INT(3, p[1].read_id);
    TEST_ASSERT_EQUAL_INT(-1, p[0].alleles[0]);
    TEST_ASSERT_EQUAL_INT(-1, p[1].alt_qi[2]);

    free_read_var_profile(p, 2);
}

static void test_update_read_var_profile_with_allele_tracks_relative_offsets(void) {
    read_var_profile_t *p = init_read_var_profile(1, 8);

    update_read_var_profile_with_allele(3, 1, 11, p);
    update_read_var_profile_with_allele(4, 0, -1, p);
    update_read_var_profile_with_allele(6, -1, -1, p);

    TEST_ASSERT_EQUAL_INT(3, p[0].start_var_idx);
    TEST_ASSERT_EQUAL_INT(6, p[0].end_var_idx);
    TEST_ASSERT_EQUAL_INT(1, p[0].alleles[0]);
    TEST_ASSERT_EQUAL_INT(11, p[0].alt_qi[0]);
    TEST_ASSERT_EQUAL_INT(0, p[0].alleles[1]);
    TEST_ASSERT_EQUAL_INT(-1, p[0].alt_qi[1]);
    TEST_ASSERT_EQUAL_INT(-1, p[0].alleles[2]);
    TEST_ASSERT_EQUAL_INT(-1, p[0].alt_qi[2]);
    TEST_ASSERT_EQUAL_INT(-1, p[0].alleles[3]);
    TEST_ASSERT_EQUAL_INT(-1, p[0].alt_qi[3]);

    free_read_var_profile(p, 1);
}

static void test_get_var_start_uses_lower_bound_binary_search(void) {
    cand_var_t vars[5] = {0};

    vars[0].pos = 10;
    vars[1].pos = 15;
    vars[2].pos = 15;
    vars[3].pos = 30;
    vars[4].pos = 40;

    TEST_ASSERT_EQUAL_INT(0, get_var_start(vars, 0, 5, 5));
    TEST_ASSERT_EQUAL_INT(1, get_var_start(vars, 0, 5, 15));
    TEST_ASSERT_EQUAL_INT(3, get_var_start(vars, 2, 5, 16));
    TEST_ASSERT_EQUAL_INT(5, get_var_start(vars, 0, 5, 41));
}

static void test_get_var_site_start_uses_lower_bound_binary_search(void) {
    var_site_t vars[5] = {0};

    vars[0].pos = 12;
    vars[1].pos = 20;
    vars[2].pos = 20;
    vars[3].pos = 21;
    vars[4].pos = 50;

    TEST_ASSERT_EQUAL_INT(0, get_var_site_start(vars, 0, 5, 1));
    TEST_ASSERT_EQUAL_INT(1, get_var_site_start(vars, 0, 5, 20));
    TEST_ASSERT_EQUAL_INT(3, get_var_site_start(vars, 2, 5, 21));
    TEST_ASSERT_EQUAL_INT(5, get_var_site_start(vars, 0, 5, 51));
}

static void test_collect_digar_from_md_tag_insertion_uses_all_inserted_base_qualities(void) {
    assert_insertion_quality_uses_all_inserted_bases(collect_digar_from_MD_tag, 1);
}

static void test_collect_digar_from_ref_seq_insertion_uses_all_inserted_base_qualities(void) {
    assert_insertion_quality_uses_all_inserted_bases(collect_digar_from_ref_seq, 0);
}

void bam_utils_suite(void) {
    RUN_TEST(test_init_read_var_profile_initializes_each_read);
    RUN_TEST(test_init_read_var_profile_with_ids_preserves_input_ids);
    RUN_TEST(test_update_read_var_profile_with_allele_tracks_relative_offsets);
    RUN_TEST(test_get_var_start_uses_lower_bound_binary_search);
    RUN_TEST(test_get_var_site_start_uses_lower_bound_binary_search);
    RUN_TEST(test_collect_digar_from_md_tag_insertion_uses_all_inserted_base_qualities);
    RUN_TEST(test_collect_digar_from_ref_seq_insertion_uses_all_inserted_base_qualities);
}
