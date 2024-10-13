#include "wavefront/wavefront_align.h"
#include "seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "abpoa.h"
#include "/homes2/yangao/programs/longcallD/edlib/include/edlib.h"


int test_edlib(char *pattern, char *text) {
    // 
    int k = -1;
    const EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    EdlibAlignResult result = edlibAlign(pattern, strlen(pattern), text, strlen(text), config);
    fprintf(stderr, "Edlib-Alignment returns score %d\n",result.editDistance);
    edlibFreeAlignResult(result);
    return 0;
}

int test_abpoa(uint8_t **bseqs, int n_seqs, int *seq_lens) {
    // Init abpoa graph
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->out_msa = 1;
    abpt->out_cons = 1;
    abpoa_post_set_para(abpt);
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, stdout);

    abpoa_cons_t *abc = ab->abc;
    for (int i = 0; i < abc->n_cons; ++i) {
        fprintf(stdout, ">Consensus_sequence");
        if (abc->n_cons > 1) {
            fprintf(stdout, "_%d ", i+1);
            for (int j = 0; j < abc->clu_n_seq[i]; ++j) { // output read ids for each cluster/group
                fprintf(stdout, "%d", abc->clu_read_ids[i][j]);
                if (j != abc->clu_n_seq[i]-1) fprintf(stdout, ",");
            }
        }
        fprintf(stdout, "\n");
        for (int j = 0; j < abc->cons_len[i]; ++j)
            fprintf(stdout, "%c", "ACGTN"[abc->cons_base[i][j]]);
        fprintf(stdout, "\n");
    }
    // Free
    abpoa_free(ab); abpoa_free_para(abpt);
    return 0;
}

int test_wfa(char *pattern, char *text) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.match = 0;
    attributes.affine_penalties.mismatch = 4;
    attributes.affine_penalties.gap_opening = 6;
    attributes.affine_penalties.gap_extension = 2;
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    // Align
    wavefront_align(wf_aligner, pattern,strlen(pattern), text,strlen(text));
    fprintf(stderr, "WFA-Alignment returns score %d\n",wf_aligner->cigar->score);
    // Display alignment
    // Free
    wavefront_aligner_delete(wf_aligner);
    return 0;
}