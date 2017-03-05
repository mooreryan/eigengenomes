#ifndef EIGG_PRINT_H
#define EIGG_PRINT_H

void
eigg_print_hashed_record(char* name,
                         char* seq,
                         char* qual,
                         unsigned long kmer_len,
                         long num_for_kmers,
                         long num_rev_kmers,
                         double* hashed_for_kmers,
                         double* hashed_rev_kmers);


void
eigg_print_hyperplanes(double** hyperplanes, long num_hyperplanes, long kmer_len);

#endif
