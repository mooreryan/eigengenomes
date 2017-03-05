#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "kseq.h"

#define VERSION "0.1"

KSEQ_INIT(gzFile, gzread)

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

void
print_hashed_record(char* name,
                    char* seq,
                    char* qual,
                    unsigned long kmer_len,
                    long num_for_kmers,
                    long num_rev_kmers,
                    double* hashed_for_kmers,
                    double* hashed_rev_kmers)
{
  int i = 0;
  printf("@%s\n", name);
  printf("%s\n+\n%s\n", seq, qual);
  printf("k, bins: [%lu", kmer_len);

  for (i = 0; i < num_for_kmers; ++i) {
    printf(",%.0f", hashed_for_kmers[i]);
  }
  for (i = 0; i < num_rev_kmers; ++i) {
    printf(",%.0f", hashed_rev_kmers[i]);
  }

  printf("]\n");
}

double
hash_encoded_kmer(double** hyperplanes,
                  int      num_hyperplanes,
                  double*  encoded_kmer,
                  int      kmer_len,
                  double*  pow2)
{
  int     hyperplane_i = 0;
  int     kmer_i = 0;
  double* this_hyperplane;
  double hash_val = 0.0;

  for (hyperplane_i = 0; hyperplane_i < num_hyperplanes; ++hyperplane_i) {
    this_hyperplane = hyperplanes[hyperplane_i];

    double sum = 0.0;

    for (kmer_i = 0; kmer_i < kmer_len - 1; ++kmer_i) {
      sum = sum + (encoded_kmer[kmer_i] * this_hyperplane[kmer_i]);
    }

    /* subtract the last...plane through origin */
    sum = sum - (encoded_kmer[kmer_i] * this_hyperplane[kmer_i]);

    if (sum > 0) {
      hash_val += pow2[hyperplane_i];
    }

    /* if (sum > 0) { */
    /*   hashed_kmer[hyperplane_i] = '1'; */
    /* } else { */
    /*   hashed_kmer[hyperplane_i] = '0'; */
    /* } */
  }

  return hash_val;
}

double* encode_ary_init()
{
  double* ary = calloc('T' + 1, sizeof(double));
  assert(ary != NULL);

  ary['A'] =  1;
  ary['C'] =  0.5;
  ary['T'] = -1;
  ary['G'] = -0.5;
  ary['N'] =  0;

  return ary;
}

void encode_ary_destroy(double* ary)
{
  free(ary);
}

double* encode_kmer(double* encoded_kmer,
                    double* encode_nt,
                    char* kmer,
                    unsigned long kmer_len)
{
  unsigned long i = 0;

  for (i = 0; i < kmer_len; ++i) {
    /* TODO assert that kmer[i] is a valid AA */
    encoded_kmer[i] = encode_nt[toupper(kmer[i])];
  }

  return encoded_kmer;
}


struct Kmers {
  char** kmers;
  long num_kmers;
  unsigned long kmer_len;
};

long get_num_kmers(long seq_len, unsigned long kmer_len)
{
  return seq_len - kmer_len + 1;
}

char** get_kmers(char* seq, unsigned long kmer_len, long num_kmers)
{
  long i = 0;
  unsigned long  j = 0;

  char** kmers = malloc(num_kmers * sizeof(char*));

  for (i = 0; i < num_kmers; ++i) {
    char* kmer = malloc((kmer_len + 1) * sizeof(char));

    for (j = 0; j < kmer_len; ++j) {
      kmer[j] = seq[i + j];
    }
    kmer[j] = '\0';

    kmers[i] = kmer;
  }

  return kmers;
}

struct Kmers* kmers_new(char* seq, long seq_len, unsigned long kmer_len)
{
  struct Kmers* kmers = malloc(sizeof(struct Kmers));
  assert(kmers != NULL);

  long num_kmers = get_num_kmers(seq_len, kmer_len);

  kmers->kmers     = get_kmers(seq, kmer_len, num_kmers);
  kmers->num_kmers = num_kmers;
  kmers->kmer_len  = kmer_len;

  return kmers;
}

void kmers_destroy(struct Kmers* kmers)
{
  long i = 0;

  for (i = 0; i < kmers->num_kmers; ++i) {
    free(kmers->kmers[i]);
  }

  free(kmers->kmers);
  free(kmers);
}

void kmers_print(struct Kmers* kmers)
{
  long i = 0;

  for (i = 0; i < kmers->num_kmers; ++i) {
    printf("#%ld -- %s\n", i, kmers->kmers[i]);
  }
}

double rand_real()
{
  /* TODO switch to random(4) for better random-ness? */
  return (double) rand() / (double) RAND_MAX;
}

/* random hyperplane through the origin */
double* hyperplane_init(int num_coef)
{
  double* hyperplane = malloc(num_coef * sizeof(double));

  int i = 0;
  double coef = 0.0;

  for (i = 0; i < num_coef; ++i) {
    coef = rand_real();

    /* if (rand_real() > 0.5) { */
    /*   coef = 1 / coef; */
    /* } */

    if (rand_real() > 0.5) {
      coef = -coef;
    }

    hyperplane[i] = coef;
  }

  return hyperplane;
}

void hyperplane_destroy(double* hyperplane)
{
  free(hyperplane);
}

/* num_coef is normally the kmer len */
double** hyperplane_ary_init(int num_hyperplanes, int num_coef)
{
  int i = 0;

  double** hyperplane_ary = malloc(num_hyperplanes * sizeof(double*));

  for (i = 0; i < num_hyperplanes; ++i) {
    hyperplane_ary[i] = hyperplane_init(num_coef);
  }

  return hyperplane_ary;
}

void hyperplane_ary_destroy(double** hyperplane_ary, int num_hyperplanes)
{
  int i = 0;

  for (i = 0; i < num_hyperplanes; ++i) {
    hyperplane_destroy(hyperplane_ary[i]);
  }

  free(hyperplane_ary);
}


void print_hyperplanes(double** hyperplanes, long num_hyperplanes, long kmer_len)
{
  for (int i = 0; i < num_hyperplanes; ++i) {
    printf("H%d: ", i);
    for (int j = 0; j < kmer_len; ++j) {
      printf("%f ", hyperplanes[i][j]);
    }
    putchar('\n');
  }
  putchar('\n');
}

/* Input a kseq, and print it out fastA/Q style */
void print_record(kseq_t* seq) {
  if (seq->qual.l) {
    putchar('@');
  } else {
    putchar('>');
  }

  printf("%s", seq->name.s);

  if (seq->comment.l) {
    printf(" %s\n", seq->comment.s);
  } else {
    putchar('\n');
  }

  printf("%s\n", seq->seq.s);

  if (seq->qual.l) { printf("+\n%s\n", seq->qual.s); }
}

int main(int argc, char* argv[])
{
  long l;
  unsigned long i = 0;
  kseq_t* seq;
  gzFile  fp;
  unsigned long num_seqs = 0;
  long seed = 0;
  unsigned long kmer_len = 0;
  unsigned long num_hyperplanes = 0;
  double** hyperplanes;
  double* encode_nt;
  double* encoded_kmer;
  double hashed_kmer;

  char* s1_name;
  char* s1_seq;
  char* s1_qual;
  unsigned long s1_seq_len = 0;

  if (argc != 5) {
    fprintf(stderr,
            "VERSION: %s\nUsage: %s "
            "<1: seed> <2: kmer size> "
            "<3: number of hyperplanes> <4: seqs.fastq>\n",
            VERSION,
            argv[0]);

    return 1;
  }

  fp = gzopen(argv[4], "r");
  if (!fp) {
    fprintf(stderr, "ERROR -- could not open %s\n", argv[4]);

    return 2;
  }

  seed = strtol(argv[1], NULL, 10);
  srand(seed);

  kmer_len = strtol(argv[2], NULL, 10);

  encoded_kmer = malloc(kmer_len * sizeof(double));

  /* Total bins = 2^num_hyperplanes */
  num_hyperplanes = strtol(argv[3], NULL, 10);

  double* pow2 = malloc(num_hyperplanes * sizeof(double));
  for (int i = 0; i < num_hyperplanes; ++i) {
    pow2[i] = pow(2, i);
  }

  seq = kseq_init(fp);
  hyperplanes = hyperplane_ary_init(num_hyperplanes, kmer_len);
  encode_nt   = encode_ary_init();

  /* print_hyperplanes(hyperplanes, num_hyperplanes, kmer_len); */

  while ((l = kseq_read(seq)) >= 0) {
    s1_name = strdup(seq->name.s);
    s1_qual = strdup(seq->qual.s);
    s1_seq  = strdup(seq->seq.s);
    s1_seq_len = seq->seq.l;

    l = kseq_read(seq);
    if (l < 0) {
      fprintf(stderr,
              "ERROR -- not enough reads in reverse file %s\n",
              argv[3]);
      return 3;
    }

    num_seqs += 2;
    if ((num_seqs % 1000) == 0) {
      fprintf(stderr, "LOG -- reading: %lu\r", num_seqs);
    }

    if (strlen(s1_seq) < kmer_len) { /* the sequence is too short for the
                                    kmer len */
      fprintf(stderr, "LOG -- Seq: %s, Len: %lu is too short\n",
              s1_name,
              s1_seq_len);

      /* do something to handle it */
      return 4;
    }

    if (seq->seq.l < kmer_len) { /* the sequence is too short for the
                                    kmer len */
      fprintf(stderr, "LOG -- Seq: %s, Len: %lu is too short\n",
              seq->name.s,
              seq->seq.l);

      /* do something to handle it */
      return 4;
    }

    struct Kmers* for_kmers = kmers_new(s1_seq, s1_seq_len, kmer_len);
    double* hashed_for_kmers = malloc(for_kmers->num_kmers * sizeof(double));

    /* the rev read might be different lenght */
    struct Kmers* rev_kmers = kmers_new(seq->seq.s, seq->seq.l, kmer_len);
    double* hashed_rev_kmers = malloc(rev_kmers->num_kmers * sizeof(double));

    for (i = 0; i < rev_kmers->num_kmers; ++i) {
      encode_kmer(encoded_kmer,
                  encode_nt,
                  rev_kmers->kmers[i],
                  kmer_len);

      hashed_kmer = hash_encoded_kmer(hyperplanes,
                                      num_hyperplanes,
                                      encoded_kmer,
                                      kmer_len,
                                      pow2);

      hashed_rev_kmers[i] = hashed_kmer;
    }

    for (i = 0; i < for_kmers->num_kmers; ++i) {
      encode_kmer(encoded_kmer,
                  encode_nt,
                  for_kmers->kmers[i],
                  kmer_len);

      hashed_kmer = hash_encoded_kmer(hyperplanes,
                                      num_hyperplanes,
                                      encoded_kmer,
                                      kmer_len,
                                      pow2);

      hashed_for_kmers[i] = hashed_kmer;
    }

    print_hashed_record(s1_name,
                        s1_seq,
                        s1_qual,
                        kmer_len,
                        for_kmers->num_kmers,
                        rev_kmers->num_kmers,
                        hashed_for_kmers,
                        hashed_rev_kmers);

    print_hashed_record(seq->name.s,
                        seq->seq.s,
                        seq->qual.s,
                        kmer_len,
                        for_kmers->num_kmers,
                        rev_kmers->num_kmers,
                        hashed_for_kmers,
                        hashed_rev_kmers);

    free(hashed_rev_kmers);
    free(hashed_for_kmers);
    free(s1_name);
    free(s1_seq);
    free(s1_qual);
    kmers_destroy(rev_kmers);
    kmers_destroy(for_kmers);
  }

  kseq_destroy(seq);
  gzclose(fp);
  hyperplane_ary_destroy(hyperplanes, num_hyperplanes);
  free(encoded_kmer);
  encode_ary_destroy(encode_nt);
  free(pow2);

  return 0;
}
