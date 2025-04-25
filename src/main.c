#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "call_var_main.h"
#include "utils.h"
#include "htslib/kstring.h"

const char PROG[20] = "longcallD";
const char DESCRIP[100] = "local-haplotagging-based small and structural variant calling";
#ifndef LONGCALLD_VERSION
const char LONGCALLD_VERSION[20] = "0.0.5";
#endif
const char CONTACT[30] = "yangao@ds.dfci.harvard.edu";
int LONGCALLD_VERBOSE = 0;
char *CMD;

static int usage(void) {//main usage
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s: %s\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", LONGCALLD_VERSION, CONTACT); 

    fprintf(stderr, "Usage:   %s <command> [options]\n\n", PROG);

    fprintf(stderr, "Command: \n");
    fprintf(stderr, "         call          call variants from long-read BAM\n");
    // fprintf(stderr, "         joint         joint variant calling for multiple samples\n");
    // fprintf(stderr, "         genotype      call genotype for given VCF\n");
    // fprintf(stderr, "         phase         phase given VCF based on long-read BAM\n");

    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[]) {
    int i; kstring_t cmd = {0,0,0};
    kputs(argv[0], &cmd);
    for (i = 1; i < argc; ++i) {
        kputs(" ", &cmd); kputs(argv[i], &cmd);
    } CMD = strdup(cmd.s); free(cmd.s);
    
    int ret = 0;
    if (argc < 2) {
        ret = 1; usage();
    } else {
        if (strcmp(argv[1], "call") == 0) ret = call_var_main(argc-1, argv+1);
        else {
            _err_error("Unrecognized command '%s'\n", argv[1]);
            ret = 1; usage();
        }
    }
    free(CMD); return ret;
}
