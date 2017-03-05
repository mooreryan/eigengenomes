#include "eigg_kmers.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

static long get_num_kmers(long seq_len, unsigned long kmer_len)
{
  return seq_len - kmer_len + 1;
}

static char** get_kmers(char* seq, unsigned long kmer_len, long num_kmers)
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

struct Kmers* eigg_kmers_new(char* seq, long seq_len, unsigned long kmer_len)
{
  struct Kmers* kmers = malloc(sizeof(struct Kmers));
  assert(kmers != NULL);

  long num_kmers = get_num_kmers(seq_len, kmer_len);

  kmers->kmers     = get_kmers(seq, kmer_len, num_kmers);
  kmers->num_kmers = num_kmers;
  kmers->kmer_len  = kmer_len;

  return kmers;
}

void eigg_kmers_destroy(struct Kmers* kmers)
{
  long i = 0;

  for (i = 0; i < kmers->num_kmers; ++i) {
    free(kmers->kmers[i]);
  }

  free(kmers->kmers);
  free(kmers);
}

void eigg_kmers_print(struct Kmers* kmers)
{
  long i = 0;

  for (i = 0; i < kmers->num_kmers; ++i) {
    printf("#%ld -- %s\n", i, kmers->kmers[i]);
  }
}
