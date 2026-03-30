#include "unity.h"

#include "bam_utils.h"

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

void bam_utils_suite(void) {
    RUN_TEST(test_init_read_var_profile_initializes_each_read);
    RUN_TEST(test_init_read_var_profile_with_ids_preserves_input_ids);
    RUN_TEST(test_update_read_var_profile_with_allele_tracks_relative_offsets);
}
