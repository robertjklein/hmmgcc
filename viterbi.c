/*
 * viterbi.h
 *
 * Module containing the actual Viterbi algorithm, as well as other algorithms
 * (forwards, backwards, etc.) as may or may not be needed.
 *
 * Robert J. Klein
 * September, 1998
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "viterbi.h"

/*
 * Returns pointer to a whole_viterbi_t structure that includes allocated
 * Viterbi and shadow arrays.  Sizes determined by n and l parameters
 * passed in (both arrays are N*L).  The arrays are unitialized.
 */
whole_viterbi_t *create_viterbi(int n_states, int l_seq) {
  whole_viterbi_t *whole_viterbi;
  int i;

  whole_viterbi = malloc(sizeof(whole_viterbi_t));  
  if (whole_viterbi == NULL) {
    fprintf (stderr, "ERROR: Could not allocate space for whole_viteribi structure\n");
    exit(2);
  }

  /* Allocate Viterbi and shadow traceback space here */
  /* Algorithm borrowed from S. Eddy in hmmer's core_algorithms.c */
  whole_viterbi->viterbi = malloc(sizeof(double)*l_seq);
  if (whole_viterbi->viterbi == NULL) {
    fprintf  (stderr, "ERROR: Could not allocate viterbi array, dimension 1\n");
    exit(2);
  }
  whole_viterbi->shadow = malloc(sizeof(int)*l_seq);
  if (whole_viterbi->shadow == NULL) {
    fprintf (stderr,"ERROR: Could not allocate shadow array, dimension 1\n");
    exit(2);
  }
  whole_viterbi->viterbi[0] = malloc(sizeof(double)*l_seq*n_states);
  if (whole_viterbi->viterbi[0] == NULL) {
    fprintf (stderr, "ERROR: Could not allocate viterbi array, dimension 2\n");
    exit(2);
  }
  whole_viterbi->shadow[0] = malloc(sizeof(int)*l_seq*n_states);
  if (whole_viterbi->shadow[0] == NULL) {
    fprintf (stderr, "ERROR: Could not allocate shadow array, dimension 2\n");
    exit(2);
  }
  for (i=1; i < l_seq; i++) {
    whole_viterbi->viterbi[i] = whole_viterbi->viterbi[0] + i * n_states;
    whole_viterbi->shadow[i] = whole_viterbi->shadow[0] + i * n_states;
  }

  /* There is no need to initalize anything as the array will be filled
     in later */
  return (whole_viterbi);
}



/*
 * Allocates and returns a "string" containing the state output of the
 * Viterbi algorithm for the given sequence.
 * Passed in the Viterbi arrays as well as the hmm parameters.
 * All probabilities are log-probabilities, as that is how they are read into
 * the hmm_t type.
 */
state_name_t *run_viterbi (whole_viterbi_t *whole_viterbi, hmm_t *hmm, char *sequence) {
  state_name_t *state_seq;
  int i_seq, i_state, i_prev_state, cur_char_as_int;
  unsigned long sequence_length;
  double max_prev;
  int state_for_max_prev;
  /* Following are for traceback */
  int final_state_max;
  double final_state_max_value;

  sequence_length = strlen(sequence);

  state_seq = malloc(sizeof(state_name_t) * sequence_length);
  if (state_seq == NULL) {
    fprintf (stderr, "ERROR: Could not allocate state sequence\n");
    exit(2);
  }

  /*
   * Initialize: For each state, the prob of being there is the sum of
   * the prob of starting there plus the appropriate emission probability
   */
  cur_char_as_int = translate_from_char(sequence[0], hmm);
  if (cur_char_as_int >= ALPHABETSIZE) {
    fprintf (stderr, "ERROR: Character %c (%d) not in alphabet %s\n", sequence[0], (int)sequence[0], hmm->alphabet);
    exit(3);
  }  
  for (i_state=0; i_state<STATES; i_state++) {
    whole_viterbi->viterbi[0][i_state] = hmm->start_probs[i_state] + hmm->emission[i_state][cur_char_as_int];
  }

  /*
   * Iterate through loop: For each state at each position, find the prob
   * of getting there from each previous state (prob from that state plus
   * transition probability).  Maximum wins, and is stored
   * in shadow array.  The prob from there is then added to the transition
   * prob and current emission prob, and put into the viterbi array.
   */
  for (i_seq = 1; i_seq < sequence_length; i_seq++) {
    cur_char_as_int = translate_from_char(sequence[i_seq], hmm);
    if (cur_char_as_int >= ALPHABETSIZE) {
      fprintf (stderr, "ERROR: Character %c (%d) not in alphabet %s\n", sequence[i_seq], (int)sequence[i_seq], hmm->alphabet);
      exit(3);
    }  
    for (i_state=0; i_state<STATES; i_state++) {
      /* First, maximize previous state */
      /* Start with state 0 as max */
      max_prev = whole_viterbi->viterbi[i_seq-1][0] + hmm->transition[0][i_state];
      state_for_max_prev = 0;
      for (i_prev_state=1; i_prev_state<STATES; i_prev_state++) {
        if (whole_viterbi->viterbi[i_seq-1][i_prev_state] + hmm->transition[i_prev_state][i_state] > max_prev) {
          max_prev = whole_viterbi->viterbi[i_seq-1][i_prev_state] + hmm->transition[i_prev_state][i_state];
          state_for_max_prev = i_prev_state;
        }
      }
      /* Now store max in the shadow array */
      whole_viterbi->shadow[i_seq][i_state] = state_for_max_prev;
      /* Add fill in the probability for this position */
      whole_viterbi->viterbi[i_seq][i_state] = max_prev + hmm->emission[i_state][cur_char_as_int];
    }
  }

  /*
   * The viterbi array is now filled in.  Now, we select the final state by
   * maximizing the P for the final state (last line of the viterbi array),
   * and use the traceback array to trace a backwards path.  Path is stored
   * in state_seq, which has already been allocated.
   */
  i_seq = sequence_length - 1;       /* Work backwards */
  /* Assume state 0 is max */
  final_state_max_value = whole_viterbi->viterbi[i_seq][0];
  final_state_max = 0;
  for (i_state = 1; i_state < STATES; i_state++) {
    if (whole_viterbi->viterbi[i_seq][i_state] > final_state_max_value) {
      final_state_max = i_state;
      final_state_max_value = whole_viterbi->viterbi[i_seq][i_state];
    }
  }
  state_seq[i_seq] = (state_name_t)final_state_max;
  while (i_seq > 0) {
    state_seq[i_seq-1] = (state_name_t)whole_viterbi->shadow[i_seq][state_seq[i_seq]];
    i_seq--;
  }
  /* Now, it's filled in.  Return i_seq */
  return(state_seq);
}

/*
 * Frees the space allocated for the Viterbi.
 * Needed for multiple-sequence FASTA files so that I don't burn up all
 * the memory needlessly.
 */
void free_viterbi (whole_viterbi_t *whole_viterbi) {
  free(whole_viterbi->viterbi[0]);
  free(whole_viterbi->shadow[0]);
  free(whole_viterbi->viterbi);
  free(whole_viterbi->shadow);
  free(whole_viterbi);
}

/*
 * This is the function that is called from hmmgcc.  This simply uses the
 * Viterbi algorithm to derive a state sequence.
 */
state_name_t *get_viterbi_state_seq (hmm_t *hmm, char *sequence) {
  state_name_t *state_seq;
  whole_viterbi_t *viterbi;

  viterbi = create_viterbi (STATES, strlen(sequence));
  state_seq = run_viterbi (viterbi, hmm, sequence);
  free_viterbi (viterbi);
  return (state_seq);
}

/* POSTERIOR DECODING */

/*
 * log_sum:
 * Given ln x and ln y, returns the sum ln (x+y).  I think a simple
 * statement of the equaltion [return (log(exp(x)+exp(y)))] gives
 * underflow errors, so I'm adapting code from sre_math.c from hmmer's
 * source tree (Specifically, DLogSum).  Specifically, I've modified it
 * to take two doubles, rather than a vector (array) and size.
 */
double log_sum (double x, double y) {
  double max = -1.0e30;
  double sum = 0.0;

  if (x > max) {
    max = x;
  }
  if (y > max) {
    max = y;
  }
  if (x > max - 50.) {
    sum += exp(x-max);
  }
  if (y > max - 50.) {
    sum += exp(y-max);
  }

  sum = log(sum) + max;
  return sum;
}

/*
 * Returns pointer to a DP matrix ready for a forward or backward run
 */
dp_matrix_t make_dp_matrix (int n_states, int l_seq) {
  dp_matrix_t dp_matrix;
  int i;

  /* Allocate DP matrix space here */
  /* Algorithm borrowed from S. Eddy in hmmer's core_algorithms.c */
  dp_matrix = malloc(sizeof(double)*l_seq);
  if (dp_matrix == NULL) {
    fprintf  (stderr, "ERROR: Could not allocate DP array, dimension 1\n");
    exit(2);
  }
  dp_matrix[0] = malloc(sizeof(double)*l_seq*n_states);
  if (dp_matrix[0] == NULL) {
    fprintf (stderr, "ERROR: Could not allocate DP array, dimension 2\n");
    exit(2);
  }
  for (i=1; i < l_seq; i++) {
    dp_matrix[i] = dp_matrix[0]+i*n_states;
  }

  /* There is no need to initalize anything as the array will be filled
     in later */
  return (dp_matrix);
}

/*
 * Do the forward or backwards algorithm.  Return value is P(x), as a log 
 * score (base e)
 */
double run_forward (dp_matrix_t forward_matrix, hmm_t *hmm, char *sequence) {
  int i_seq, i_state, i_prev_state, cur_char_as_int;
  unsigned long sequence_length;
  double sum_over_states;

  sequence_length = strlen(sequence);
  /*
   * Initialize: For each state, the prob of being there is the sum of
   * the probaility of starting there plus the appropriate emission probability
   */
  cur_char_as_int = translate_from_char(sequence[0], hmm);
  if (cur_char_as_int >= ALPHABETSIZE) {
    exit(3);
  }  
  for (i_state=0; i_state<STATES; i_state++) {
    forward_matrix[0][i_state] = hmm->start_probs[i_state] + hmm->emission[i_state][cur_char_as_int];
  }

  /*
   * Iterate through loop: For each state at each position, find the prob
   * of getting there from each previous state (prob from that state plus
   * transition probability).  The sum of this is taken, and added to the
   * emission prob.  This is then stored in the forward array
   */
  for (i_seq = 1; i_seq < sequence_length; i_seq++) {
    cur_char_as_int = translate_from_char(sequence[i_seq], hmm);
    if (cur_char_as_int >= ALPHABETSIZE) {
      fprintf (stderr, "ERROR: Character %c (%d) not in alphabet %s\n", sequence[i_seq], (int)sequence[i_seq], hmm->alphabet);
      exit(3);
    }  
    for (i_state=0; i_state<STATES; i_state++) {
      /* First, set the total to the result from state 0 */
      sum_over_states = forward_matrix[i_seq-1][0] + hmm->transition[0][i_state];
      for (i_prev_state=1; i_prev_state<STATES; i_prev_state++) {
	sum_over_states = log_sum (sum_over_states, (forward_matrix[i_seq-1][i_prev_state] + hmm->transition[i_prev_state][i_state]));
      }
      /* Now add the emission prob, and store in forward matrix */
      forward_matrix[i_seq][i_state] = sum_over_states + hmm->emission[i_state][cur_char_as_int];
    }
  }

  /* The matrix is now filled in, so calculate P(x) and return it */
  /* Assumes a sub k0 (probability of going from state to null state) is 1/L */
  sum_over_states = forward_matrix[sequence_length-1][0] + log(1.0/(double)sequence_length);
  for (i_state=1; i_state<STATES; i_state++) {
    sum_over_states = log_sum (sum_over_states, forward_matrix[sequence_length-1][i_state] + log(1.0/(double)sequence_length));
  }
  return (sum_over_states);
}

/*
 * Same code as forward, with appropriate modifications
 */
double run_backward (dp_matrix_t backward_matrix, hmm_t *hmm, char *sequence) {
  int i_seq, i_state, i_next_state, cur_char_as_int;
  unsigned long sequence_length;
  double sum_over_states;

  sequence_length = strlen(sequence);

  /*
   * Initialize: For each end state, the transition probability between there
   * and state "0" (a null, start state) is 1/sequnece_length
   */
  for (i_state = 0; i_state<STATES; i_state++) {
    backward_matrix[sequence_length-1][i_state] = log (1.0/(double)sequence_length);
  }

  /*
   * Iterate through loop: For each state for the "next" letter (i+1), take
   * the sum of the log probabilities of the transition, the emission for
   * the "next" letter, and the backwards value for the "next" letter
   */
  for (i_seq = sequence_length-2; i_seq >=0; i_seq--) {
    cur_char_as_int = translate_from_char(sequence[i_seq+1], hmm);
    if (cur_char_as_int >= ALPHABETSIZE) {
      fprintf (stderr, "ERROR: Character %c (%d) not in alphabet %s\n", sequence[i_seq], (int)sequence[i_seq], hmm->alphabet);
      exit(3);
    } 
    for (i_state=0; i_state<STATES; i_state++) {
      /* First, set the total to the result from state 0 */
      sum_over_states = backward_matrix[i_seq+1][0] + hmm->transition[i_state][0] + hmm->emission[0][cur_char_as_int];
      for (i_next_state=1; i_next_state<STATES; i_next_state++) {
	sum_over_states = log_sum (sum_over_states, (backward_matrix[i_seq+1][i_next_state] + hmm->transition[i_state][i_next_state] + hmm->emission[i_next_state][cur_char_as_int]));
      }
      /* Now store in backward matrix */
      backward_matrix[i_seq][i_state] = sum_over_states;
    }
  }

  /* The matrix is now filled in, so calculate P(x) and return it */
  cur_char_as_int = translate_from_char(sequence[0], hmm);
  if (cur_char_as_int >= ALPHABETSIZE) {
    fprintf (stderr, "ERROR: Character %c (%d) not in alphabet %s\n", sequence[0], (int)sequence[0], hmm->alphabet);
    exit(3);
  } 
  sum_over_states = backward_matrix[0][0] + hmm->start_probs[0] + hmm->emission[0][cur_char_as_int];
  for (i_state=1; i_state<STATES; i_state++) {
    sum_over_states = log_sum (sum_over_states, backward_matrix[0][i_state] + hmm->start_probs[i_state] + hmm->emission[i_state][cur_char_as_int]);
  }
  return (sum_over_states);
}


/* Returns the posterior decoding for the given state and parameters */
double *posterior_decode (dp_matrix_t forward_matrix, dp_matrix_t backward_matrix, double P, unsigned long sequence_length, state_name_t state) {
  double *prob_seq;
  int i;

  prob_seq = malloc(sizeof(double) * sequence_length);
  if (prob_seq == NULL) {
    fprintf (stderr, "ERROR: Could not allocate posterior probability array\n");
    exit(2);
  }

  for (i = 0; i<sequence_length; i++) {
    prob_seq[i] = forward_matrix[i][state] + backward_matrix[i][state] - P;
  }

  return (prob_seq);
}

/*
 * Given a cutoff, returns a state seq for it.  Called both from hmmgcc
 * (for -p) as well as from below function (for -P)
 */
state_name_t *get_posterior_cutoff_state_seq (hmm_t *hmm, char *sequence, double cutoff, state_name_t over_state, state_name_t under_state) {
  dp_matrix_t forward, backward;
  double Pf, Pb, P;
  double *decoded;
  state_name_t *state_seq;
  unsigned long sequence_length;
  unsigned long c;

  sequence_length = strlen(sequence);

  /* First, do forward and backwards algorithms */
  forward = make_dp_matrix (STATES, sequence_length);
  backward = make_dp_matrix (STATES, sequence_length);
  Pf = run_forward (forward, hmm, sequence);
  Pb = run_backward (backward, hmm, sequence);
  P = (Pf + Pb) / 2.0;

  /* Now, do posterior decoding */
  decoded = posterior_decode (forward, backward, P, sequence_length, over_state);

  /* And create a state sequence from it */
  state_seq = malloc(sizeof(state_name_t) * sequence_length);
  if (state_seq == NULL) {
    fprintf (stderr, "ERROR: Could not allocate state sequence\n");
    exit(2);
  }

  for (c=0; c<sequence_length; c++) {
    if (decoded[c] >= cutoff) {
      state_seq[c] = over_state;
    } else {
      state_seq[c] = under_state;
    }
  }

  /* And free everything */
  free(forward[0]);
  free(backward[0]);
  free(forward);
  free(backward);
  free(decoded);

  return (state_seq);
}

/*
 * Second method.  Given an array with ranges, gets state seq with cutoff that
 * insures everything in the ranges is included.
 */
state_name_t *get_posterior_greedy_state_seq (hmm_t *hmm, char *sequence, int *include, state_name_t true_state, state_name_t false_state) {
  /* Based on code above for posterior_cutoff */
  dp_matrix_t forward, backward;
  double Pf, Pb, P;
  double *decoded;
  state_name_t *state_seq;
  unsigned long sequence_length;
  unsigned long c;
  double cutoff;

  sequence_length = strlen(sequence);

  /* First, do forward and backwards algorithms */
  forward = make_dp_matrix (STATES, sequence_length);
  backward = make_dp_matrix (STATES, sequence_length);
  Pf = run_forward (forward, hmm, sequence);
  Pb = run_backward (backward, hmm, sequence);
  P = (Pf + Pb) / 2.0;

  /* Now, do posterior decoding */
  decoded = posterior_decode (forward, backward, P, sequence_length, true_state);

  /* Here is the intermediate step: determine a cutoff such that every position
     in include comes in above the cutoff */
  cutoff = 1;
  for (c=0; c<sequence_length; c++) {
    if (include[c] == 1 && decoded[c] < cutoff) {
      cutoff = decoded[c];
    }
  }

  /* Print out the cutoff to stderr for reference purposes */
  fprintf (stderr, "Cutoff = %f, ln(Cutoff) = %f\n", exp(cutoff), cutoff);

  /* And create a state sequence from it */
  state_seq = malloc(sizeof(state_name_t) * sequence_length);
  if (state_seq == NULL) {
    fprintf (stderr, "ERROR: Could not allocate state sequence\n");
    exit(2);
  }

  for (c=0; c<sequence_length; c++) {
    if (decoded[c] >= cutoff) {
      state_seq[c] = true_state;
    } else {
      state_seq[c] = false_state;
    }
  }

  /* And free everything */
  free(forward[0]);
  free(backward[0]);
  free(forward);
  free(backward);
  free(decoded);

  return (state_seq);
}


/*
 * Third method.  Simply print out the posterior decoding log probabilities
 * along with the sequence position and sequence member.
 * Doesn't return anything.
 */
void print_all_posterior (hmm_t *hmm, char *sequence, state_name_t true_state, char *seqname) {
  /* Based on code above for posterior_cutoff */
  dp_matrix_t forward, backward;
  double Pf, Pb, P;
  double *decoded;
  state_name_t *state_seq;
  unsigned long sequence_length;
  unsigned long c;

  sequence_length = strlen(sequence);

  /* First, do forward and backwards algorithms */
  forward = make_dp_matrix (STATES, sequence_length);
  backward = make_dp_matrix (STATES, sequence_length);
  Pf = run_forward (forward, hmm, sequence);
  Pb = run_backward (backward, hmm, sequence);
  P = (Pf + Pb) / 2.0;

  /* Now, do posterior decoding */
  decoded = posterior_decode (forward, backward, P, sequence_length, true_state);

  /* Interate through, and print results */
  for (c=0; c<sequence_length; c++) {
    printf ("%s\t%d\t%c\t%.2f\n", seqname, c+1, sequence[c], decoded[c]);
  }

  /* And free everything */
  free(forward[0]);
  free(backward[0]);
  free(forward);
  free(backward);
  free(decoded);

}

