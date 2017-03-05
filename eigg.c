#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "kseq.h"
#include "eigg_encode.h"
#include "eigg_hash.h"
#include "eigg_hyperplane.h"
#include "eigg_kmers.h"
#include "eigg_print.h"

#define VERSION "0.1"

KSEQ_INIT(gzFile, gzread)

/* Input a kseq, and print it out fastA/Q style */
void print_record(kseq_t* seq)
{
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
  hyperplanes = eigg_hyperplane_ary_init(num_hyperplanes, kmer_len);
  encode_nt   = eigg_encode_ary_init();

  /* eigg_print_hyperplanes(hyperplanes, num_hyperplanes, kmer_len); */

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

    struct Kmers* for_kmers = eigg_kmers_new(s1_seq, s1_seq_len, kmer_len);
    double* hashed_for_kmers = malloc(for_kmers->num_kmers * sizeof(double));

    /* the rev read might be different lenght */
    struct Kmers* rev_kmers = eigg_kmers_new(seq->seq.s, seq->seq.l, kmer_len);
    double* hashed_rev_kmers = malloc(rev_kmers->num_kmers * sizeof(double));

    for (i = 0; i < rev_kmers->num_kmers; ++i) {
      eigg_encode_kmer(encoded_kmer,
                       encode_nt,
                       rev_kmers->kmers[i],
                       kmer_len);

      hashed_kmer = eigg_hash_encoded_kmer(hyperplanes,
                                           num_hyperplanes,
                                           encoded_kmer,
                                           kmer_len,
                                           pow2);

      hashed_rev_kmers[i] = hashed_kmer;
    }

    for (i = 0; i < for_kmers->num_kmers; ++i) {
      eigg_encode_kmer(encoded_kmer,
                       encode_nt,
                       for_kmers->kmers[i],
                       kmer_len);

      hashed_kmer = eigg_hash_encoded_kmer(hyperplanes,
                                           num_hyperplanes,
                                           encoded_kmer,
                                           kmer_len,
                                           pow2);

      hashed_for_kmers[i] = hashed_kmer;
    }

    eigg_print_hashed_record(s1_name,
                             s1_seq,
                             s1_qual,
                             kmer_len,
                             for_kmers->num_kmers,
                             rev_kmers->num_kmers,
                             hashed_for_kmers,
                             hashed_rev_kmers);

    eigg_print_hashed_record(seq->name.s,
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
    eigg_kmers_destroy(rev_kmers);
    eigg_kmers_destroy(for_kmers);
  }

  kseq_destroy(seq);
  gzclose(fp);
  eigg_hyperplane_ary_destroy(hyperplanes, num_hyperplanes);
  free(encoded_kmer);
  eigg_encode_ary_destroy(encode_nt);
  free(pow2);

  return 0;
}
