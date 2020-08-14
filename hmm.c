/*
 * hmm.c
 *
 * Functions and routines for general manipulation of the HMM.
 * Includes initialization routines for the HMM, so HMM definition can
 * get modified here as well.
 *
 * Robert J. Klein
 * September, 1998
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hmm.h"

/*
 * This function returns a character randomly chosen from the string.
 * Used in conjunction with the function below that resolves degenerate code
 * nucleotides.
 */
char random_from_string (char *s) {
  int i;
  do {
    i = (int) ((float)(strlen(s)-1)*rand()/(RAND_MAX+1.0));
  } while (i<0 || i>=strlen(s));
  return(s[i]);
}

/*
 * This function resolves "degnerate" nucleotides by selecting a random 
 * A, C, G, or T as appropriate by the code present there.  Returns
 * the character passed in if that character does not represent a
 * non-degnerate nucleotide (either A, C, G, or T or not representative
 * at all of a nucleotide.
 *
 * The degenerate code used here is:
 * (taken from http://www.neb.com/neb/products/REs/RE_code.html
 *
 *                         R = G or A
 *                         K = G or T
 *                         B = not A (C or G or T)
 *                         V = not T (A or C or G)
 *                         Y = C or T
 *                         S = G or C
 *                         D = not C (A or G or T)
 *                         N = A or C or G or T
 *                         M = A or C
 *                         W = A or T
 *                         H = not G (A or C or T)
 *
 * This function assumes all letters are already uppercased via toupper
 * before calling.  In other words, it will return a "n" if passed an "n"
 * because it will assume that the symbol for all nucleotides will be passed
 * in as "N".
 */
char resolve_degenerate (char c) {
  c = toupper(c);
  switch (c) {
    case 'R' : return(random_from_string("GA"));
    case 'K' : return(random_from_string("GT"));
    case 'B' : return(random_from_string("CGT"));
    case 'V' : return(random_from_string("ACG"));
    case 'Y' : return(random_from_string("CT"));
    case 'S' : return(random_from_string("GC"));
    case 'D' : return(random_from_string("AGT"));
    case 'N' : return(random_from_string("ACGT"));
    case 'M' : return(random_from_string("AC"));
    case 'W' : return(random_from_string("AT"));
    case 'H' : return(random_from_string("ACT"));
  }
  return(c);
}

/*
 * Translate from the char to the position it is in the alphabet array
 * Returning an out-of-bounds number means the character isn't in the
 * alphabet.  First substitutes an appropriate random character for characters
 * representing degenerate positions.  That is done via the function above.
 */
int translate_from_char (char c, hmm_t *hmm) {
  int i = 0;
  c = resolve_degenerate(toupper(c));
  while (hmm->alphabet[i] != '\0' && hmm->alphabet[i] != c) i++;
  return(i);
}

/* 
 * Returns a pointer to a newly allocated hmm that is filled with emission
 * and transition probabilities.  Probabilities are read from the file
 * filename.
 */
hmm_t *create_hmm (char *filename) {
  hmm_t *hmm;
  int i,j,k;
  FILE *f;
  char buffer[256];
  char *ch;
  /* 
   * This array stores the doubles as read in from the HMM file, as they
   * are converted to logs.  The array is then scanned down, with the data
   * filled in as appropriate.
   */  
  double temp_array[STATES*ALPHABETSIZE+STATES*STATES+2];
  

  hmm = malloc(sizeof(hmm_t));
  if (hmm == NULL) {
    fprintf (stderr, "ERROR: Couldn't allocate space for hmm_t\n");
    exit(1);
  }

  /* Allocate space for the alphabet string */
  hmm->alphabet = malloc(strlen(ALPHABET)+1);
  strcpy (hmm->alphabet, ALPHABET);

  /* Allocate emission probability space here */
  /* Algorithm borrowed from S. Eddy in hmmer's core_algorithms.c */
  hmm->emission = malloc(sizeof(double)*STATES);
  hmm->emission[0] = malloc(sizeof(double)*STATES*ALPHABETSIZE);
  for (i=1; i < STATES; i++) {
    hmm->emission[i] = hmm->emission[0] + i * ALPHABETSIZE;
  }

  /* Allocate transition probability space here */
  hmm->transition = malloc(sizeof(double)*STATES);
  hmm->transition[0] = malloc(sizeof(double)*STATES*STATES);
  for (i=1; i < STATES; i++) {
    hmm->transition[i] = hmm->transition[0] + i * STATES;
  }

  /* Read in the data */
  f = fopen(filename,"r");
  if (f == NULL) {
    fprintf (stderr, "ERROR: Could not open hmm data file %s\n", filename);
    exit(2);
  }
  i = 0;
  while (fgets(buffer, 255, f)) {
    if (buffer[0] == '#') {
      continue;
    }
    temp_array[i++] = log(atof(buffer)); 
  }
  k=0;
  hmm->start_probs[GC] = temp_array[k++];
  hmm->start_probs[AT] = temp_array[k++];
  for (i=0; i<STATES; i++) 
    for (j=0; j<ALPHABETSIZE; j++) 
      hmm->emission[i][j] = temp_array[k++];
  for (i=0; i<STATES; i++)
    for (j=0; j<STATES; j++)
      hmm->transition[i][j] = temp_array[k++];
  (void)fclose(f);

  return(hmm);
}
  





