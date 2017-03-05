#include "eigg_print.h"
#include <stdio.h>

void
eigg_print_hashed_record(char* name,
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

void eigg_print_hyperplanes(double** hyperplanes, long num_hyperplanes, long kmer_len)
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
