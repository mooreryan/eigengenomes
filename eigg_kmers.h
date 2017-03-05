#ifndef EIGG_KMERS_H
#define EIGG_KMERS_H

struct Kmers {
  char** kmers;
  long num_kmers;
  unsigned long kmer_len;
};

struct Kmers* eigg_kmers_new(char* seq, long seq_len, unsigned long kmer_len);

void eigg_kmers_destroy(struct Kmers* kmers);

void eigg_kmers_print(struct Kmers* kmers);

#endif
