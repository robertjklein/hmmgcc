/*
 * hmmgcc.c
 *
 * Main program for finding GC content using Hidden Markov Models.  Passed
 * in two parameters -- the filename of the sequence and the filename of the
 * HMM -- and returns the optimal (Viterbi) state sequence on stdout.  Other
 * command-line options allow alternative hmm algorithms (e.g. forward-
 * backward).
 *
 * Robert J. Klein
 * September, 1998
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "hmm.h"
#include "viterbi.h"

typedef enum { unknown, viterbi_mode, posterior_cutoff_mode, posterior_greedy_mode, posterior_all_mode } hmmgcc_mode_t;

/*
 * Prints out usage message on stderr and exits
 */
void usage () {
  fprintf (stderr, "USAGE: hmmgcc [ -l length ] [ -v | -p cutoff | -P range-file | -a ] seq-file hmm-file\n");
  fprintf (stderr, "\t-v\tViterbi mode -- print Viterbi state sequence (DEFAULT)\n");
  fprintf (stderr, "\t-p\tPosterior cutoff mode -- uses provided cutoff (ln P)\n");
  fprintf (stderr, "\t-P\tPosterior \"greedy\" mode -- uses highest possible cutoff to\n");
  fprintf (stderr, "\t\tinclude all areas specified in the range file\n");
  fprintf (stderr, "\t-a\tPosterior \"print all\" mode -- prints decoding for all nts\n");
  fprintf (stderr, "\t-l\tAll members of state looking for must be at least l long\n");
  fprintf (stderr, "\n\n");
  exit (1);
}



/*
 * Reads in the next fasta sequence in FILE *.  Allocation done on the pointers
 * passed in by reference (assigned at end).  Returns 0 if no sequence left 
 * in file, 1 if there is a sequence that was read in.
 */
int read_next_sequence (FILE *f, char **name_return, char **seq_return) {
  int cur_allocated_size;
  char *fullseq, *name;
  int pos;
  int ch;

  /* Start by allowing for 4MB of sequence */
  cur_allocated_size = 4194304;
  fullseq = malloc(cur_allocated_size * sizeof(char));
  if (fullseq == NULL) {
    fprintf (stderr, "ERROR: Could not allocate memory for sequence\n");
    exit(1);
  }

  /* Allocated 40 characters max for the name of the seq */
  name = malloc(40);
  if (name == NULL) {
    fprintf (stderr, "ERROR: Could not allocate memory for name \n");
    exit(1);
  }

  /* Read in the name */
  if (fgets (name, 40-1, f) == NULL ) {
    return (0);                 /* Not a complete sequence */
  }
  if (strchr (name, (int)'\n') == NULL) {     /* No newline here */
    do {
      ch = fgetc(f);
    } while (ch != '\n');
  }
  for (pos=0; !isspace(name[pos]) && pos < 40; pos++);
  name[pos] = '\0';
  /* 
   * Thus, the name is the shorter of the string before the first space or
   * the first 39 characters
   */

  /* Now, scan in the sequence */
  pos = 0;
  do {
    ch = fgetc(f);
    if (ch == '\n' || ch == '-') continue;	
      /* Skip "gaps" -- shouldn't ever be there, but exists at least 
      in S. cerevisiae FASTA files, and probably others as well */
    if (ch == EOF || ch == '>') break;
    if (cur_allocated_size <= pos) {
      cur_allocated_size += 1048576;        /* Increase the size by 1MB */
      fullseq = realloc (fullseq, cur_allocated_size);
      if (fullseq == NULL) {
	fprintf (stderr, "ERROR: Could not fully allocated sequence\n");
	exit(1);
      }
    }
    fullseq[pos] = (char)ch;
    pos++;
  } while (1);
  if (cur_allocated_size <= pos) {
    cur_allocated_size += 1;
    fullseq = realloc (fullseq, cur_allocated_size);
    if (fullseq == NULL) {
      fprintf (stderr, "ERROR: Could not fully allocate sequence by 1 byte\n");
      exit(1);
    }
  }
  fullseq[pos] = '\0';
  /* If memory becomes tight, do a realloc here so that no more than absolutely
   * necessary is used */

  *name_return = name;
  *seq_return = fullseq;
  if (fullseq[0] == '\0') {
    return(0);
  } else {
    return(1);
  }
}

/*
 * Creates the greedy list, and marks as "True" those ranges listed in
 * the greedy file
 */
int *make_greedy_list (char *greedy_file, char *sequence_name, unsigned long sequence_length) {
  unsigned long c;
  FILE *f;
  int *list;
  char buffer[1024];
  char *cp;
  unsigned long start, finish;

  list = malloc(sizeof(int)*sequence_length);
  if (list == NULL ) {
    fprintf (stderr, "ERROR: Cannot allocate space for greedy list\n");
    exit(2);
  }
  for (c=0; c<sequence_length; c++) {
    list[c] = 0;
  }

  f = fopen (greedy_file, "r");
  if (f == NULL) {
    fprintf (stderr, "ERROR: Could not open file %s\n", greedy_file);
    exit(2);
  }

  while (fgets(buffer, 1023, f)) {
    cp = buffer;
    while (!isspace(*cp)) cp++;        /* Until first whitespace */
    *cp = '\0';                        /* Set to NULL for later */
    cp++;
    while (isspace(*cp)) cp++;         /* Where the first number is */
    start = atoi(cp) - 1;
    while (*cp != '-') cp++;
    cp++;
    finish=atoi(cp) - 1;
    if (start > finish) {
      c = start;
      start = finish;
      finish = c;
    }
    if (strcmp (buffer, sequence_name) == 0) {
      for (c=start; c<=finish; c++) {
	list[c] = 1;
      }
    }
  }
  (void) fclose(f);

  return (list);
}
  

/*
 * Print out the state sequence passed in.  Only prints the state represented
 * by the true_state variable, and only if length is at least min_length.
 * Do it elegantly, by listing ranges of positions for each state.  Add 1 to 
 * the array indices before printing so that actual sequence position as a 
 * biologist sees it will be displayed.  
 */
void print_state_seq (state_name_t *seq, int len, char *sequence_name, state_name_t true_state, int min_len) {
  int i;
  int cur_start;
  state_name_t cur_state;

  cur_start = 0;
  cur_state = seq[0];

  for (i=1; i<len; i++) {
    if (seq[i] != cur_state) {
      if (cur_state == true_state && (i-(cur_start+1) + 1) >= min_len) {
        printf ("%s\t%d-%d\n", sequence_name, cur_start+1, i);
      }
      cur_state = seq[i];
      cur_start = i;
    }
  }
  if (cur_state == true_state && (len-(cur_start+1) + 1) >= min_len) {
    printf ("%s\t%d-%d\n", sequence_name, cur_start+1, len);
  }
}

/*
 * Two arguments -- first is the sequence filename, the second is the HMM
 * filename
 */
int main (int argc, char *argv[]) {
  char *sequence, *sequence_name;
  int c, ch;
  double cutoff;
  hmm_t *hmm;
  state_name_t *state_seq;
  FILE *sequence_file;
  hmmgcc_mode_t mode = unknown;
  int *greedy_list;
  char *greedy_filename = NULL;
  char *hmm_filename = NULL;
  char *sequence_filename = NULL;
  int min_length = 0;               /* Default is no minimum */

  for (c=1; c<argc; c++) {
    if (argv[c][0] == '-') {
      switch (argv[c][1]) {
      case 'v': 
	mode = viterbi_mode;
	break;
      case 'p':
	mode = posterior_cutoff_mode;
	c++;
	cutoff = atof (argv[c]);
	break;
      case 'P':
	mode = posterior_greedy_mode;
	c++;
	greedy_filename = argv[c];
	break;
      case 'a':
	mode = posterior_all_mode;
	break;
      case 'l':
	c++;
	min_length = atof(argv[c]);
	break;
      default:
	usage();
      }
    } else if (sequence_filename == NULL) {
      sequence_filename = argv[c];
    } else if (hmm_filename == NULL) {
      hmm_filename = argv[c];
    } else {
      usage();
    }
  }

  /* If mode not set, use default */
  if (mode == unknown) mode=viterbi_mode;
  if (hmm_filename == NULL || sequence_filename == NULL) usage();

  sequence_file = fopen (sequence_filename, "r");
  if (sequence_file == NULL) {
    fprintf (stderr, "ERROR: Could not open %s\n", sequence_filename);
    exit(1);
  }
  do {
    ch = fgetc(sequence_file);
  } while (ch != '>' && ch != EOF);

  hmm = create_hmm(hmm_filename);

  while (read_next_sequence (sequence_file, &sequence_name, &sequence)) {
    switch (mode) {
    case viterbi_mode:
      state_seq = get_viterbi_state_seq (hmm, sequence);
      break;
    case posterior_cutoff_mode:
      state_seq = get_posterior_cutoff_state_seq (hmm, sequence, cutoff, GC, AT);
      break;
    case posterior_greedy_mode:
      /* First read in the ranges for the current sequence_name */
      greedy_list = make_greedy_list (greedy_filename, sequence_name, strlen(sequence));
      state_seq = get_posterior_greedy_state_seq (hmm, sequence, greedy_list, GC, AT);
      free(greedy_list);
      break;
    case posterior_all_mode:
      print_all_posterior (hmm, sequence, GC, sequence_name);
      break;
    }
    /* 
     * Prints out the state sequence.  This can be modified to be as nice or
     * as messy as I want. 
     */
    if (mode != posterior_all_mode) {
      print_state_seq (state_seq, strlen(sequence), sequence_name, GC, min_length);
      free(state_seq);
    }
    free(sequence_name);
    free(sequence);
  }

  (void)fclose(sequence_file);

}

  














