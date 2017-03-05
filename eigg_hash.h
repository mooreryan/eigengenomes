#ifndef EIGG_HASH_H
#define EIGG_HASH_H

double
eigg_hash_encoded_kmer(double** hyperplanes,
                       int      num_hyperplanes,
                       double*  encoded_kmer,
                       int      kmer_len,
                       double*  pow2);

#endif
