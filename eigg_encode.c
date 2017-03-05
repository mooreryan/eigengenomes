#include "eigg_encode.h"

#include <assert.h>
#include <stdlib.h>
#include <ctype.h>

double* eigg_encode_ary_init()
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

void eigg_encode_ary_destroy(double* ary)
{
  free(ary);
}

double* eigg_encode_kmer(double* encoded_kmer,
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
