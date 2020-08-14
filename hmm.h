/*
 * hmm.h
 *
 * Data structures for use with Hidden Markov models.
 *
 * This file contains several numbers that get modified as the model changes.
 * (Specifically, alphabet and STATES)
 *
 * Robert J. Klein
 * September, 1998
 */

#ifndef _hmm_h
#define _hmm_h

/* Define alphabet */
/* Note that this isn't exactly precise, as hmm.c contains functions
 * for randomly choosing a base when it's a degenerate code in the
 * sequence data */
#define ALPHABETSIZE 4
#define ALPHABET "ACGT"

/* Define states */
#define STATES 2
typedef enum { GC, AT } state_name_t;

/* 
 * Probability arrays are dynamically allocated via malloc.  This just
 * defines them for convenience
 */
typedef double **emission_t;
typedef double **transition_t;

typedef struct _hmm_t {
  emission_t emission;
  transition_t transition;
  double start_probs[STATES];
  char *alphabet;
} hmm_t;


/*
 * FUNCTIONS BELOW
 */

/* Translate from char to int for the alphabet */
int translate_from_char (char c, hmm_t *hmm);

/* Allocates space for, and fills probabilities, of an hmm */
hmm_t *create_hmm (char *filename);

#endif


