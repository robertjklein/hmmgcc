/*
 * viterbi.h
 *
 * Module containing the actual Viterbi algorithm, as well as other algorithms
 * (forwards, backwards, etc.) as may or may not be needed.
 *
 * Robert J. Klein
 * September, 1998
 */

#ifndef _viterbi_h
#define _viterbi_h

#include "hmm.h"

/* DP Matrix and traceback arrays */
typedef double **dp_matrix_t;
typedef int **shadow_t;
typedef struct _whole_viterbi_t {
  dp_matrix_t viterbi;
  shadow_t shadow;
} whole_viterbi_t;

/*
 * Returns the state sequence for a given sequence and HMM using the
 * Viterbi algorithm.
 */
state_name_t *get_viterbi_state_seq (hmm_t *hmm, char *sequence);

/* POSTERIOR DECODING */

/*
 * First method.  Given a cutoff, returns a state seq for it
 */
state_name_t *get_posterior_cutoff_state_seq (hmm_t *hmm, char *sequence, double cutoff, state_name_t over_state, state_name_t under_state);

/*
 * Second method.  Given an array with ranges, gets state seq with cutoff that
 * insures everything in the ranges is included.
 */
state_name_t *get_posterior_greedy_state_seq (hmm_t *hmm, char *sequence, int *include, state_name_t true_state, state_name_t false_state);

/*
 * Third method.  Simply dumps the posterior decoding to STDOUT
 */
void print_all_posterior (hmm_t *hmm, char *sequence, state_name_t true_state, char *seqname);

#endif




