#include "unity.h"

const char PROG[20] = "longcallD";
int LONGCALLD_VERBOSE = 0;
char *CMD = (char *)"test/unit/unit_tests";

void bam_utils_suite(void);
void collect_var_suite(void);

void setUp(void) {}
void tearDown(void) {}

int main(void) {
    UNITY_BEGIN();
    bam_utils_suite();
    collect_var_suite();
    return UNITY_END();
}
