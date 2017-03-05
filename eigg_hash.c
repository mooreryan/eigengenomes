#include "eigg_hash.h"

double
eigg_hash_encoded_kmer(double** hyperplanes,
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
