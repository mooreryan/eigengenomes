#ifndef EIGG_ENCODE_H
#define EIGG_ENCODE_H

double* eigg_encode_ary_init();

void eigg_encode_ary_destroy(double* ary);

double* eigg_encode_kmer(double* encoded_kmer,
                         double* encode_nt,
                         char* kmer,
                         unsigned long kmer_len);

#endif
