#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "call_var.h"
#include "htslib/kstring.h"
// #include "collect_snps.h"

const char PROG[20] = "longcallD";
const char DESCRIP[100] = "";
const char VERSION[20] = "0.0.1";
const char CONTACT[30] = "gaoy1@chop.edu";
int LONGCALLD_VERBOSE = 0;
char *CMD;

static int usage(void)	//main usage
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: %s%s\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", VERSION, CONTACT); 

	fprintf(stderr, "Usage:   %s <command> [options]\n\n", PROG);

	fprintf(stderr, "Command: \n");
    fprintf(stderr, "         call          call variants from long-read BAM\n");
	// fprintf(stderr, "         index              index reference graph based on GFA file\n");

	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
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
		// else if (strcmp(argv[1], "index") == 0)   return index_main(argc, argv);
		// else if (strcmp(argv[1], "map") == 0)   return map_main(argc, argv);
		else {
			fprintf(stderr, "[main] Unrecognized command '%s'\n", argv[1]);
			ret = 1; usage();
		}
	}
	free(CMD); return ret;
}