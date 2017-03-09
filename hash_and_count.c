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
#include "eigg_version.h"

#include "vendor/tommyarray.h"
#include "vendor/tommyhashlin.h"

KSEQ_INIT(gzFile, gzread)

struct hash_bucket_count_t {
  unsigned long hash_bucket;
  unsigned long count;

  tommy_node node;
};

int hash_bucket_count_compare(const void* arg, const void* hash_bucket_count)
{
  return *(const unsigned long*)arg !=
    ((const struct hash_bucket_count_t*)hash_bucket_count)->hash_bucket;
}

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
  unsigned long kmer_len = 0;
  unsigned long num_seen_hash_buckets = 0;
  unsigned long num_hyperplanes = 0;
  unsigned long num_seqs = 0;
  long seed = 0;

  kseq_t* seq;
  gzFile  fp;

  double** hyperplanes;
  double* encode_nt;
  double* encoded_kmer;

  unsigned long hashed_kmer;

  char* s1_name;
  char* s1_seq;
  char* s1_qual;
  unsigned long s1_seq_len = 0;

  tommy_array*   hash_bucket_names = malloc(sizeof(tommy_array));
  assert(hash_bucket_names != NULL);
  tommy_hashlin* hash_bucket_counts = malloc(sizeof(tommy_hashlin));
  assert(hash_bucket_counts != NULL);

  /* to hold the queries of the hash table */
  struct hash_bucket_count_t* hash_bucket_count;

  if (argc != 5) {
    fprintf(stderr,
            "VERSION: %s\nUsage: %s "
            "<1: seed> <2: kmer size> "
            "<3: number of hyperplanes> <4: seqs.fastq> "
            "> seqs.hash_counts\n",
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

  /* TODO toss out the first few random values */

  kmer_len = strtol(argv[2], NULL, 10);

  encoded_kmer = malloc(kmer_len * sizeof(double));
  assert(encoded_kmer != NULL);

  /* Total bins = 2^num_hyperplanes */
  num_hyperplanes = strtol(argv[3], NULL, 10);


  unsigned long num_hash_buckets =
    (unsigned long)pow(2, num_hyperplanes);

  /* unsigned long* hashed_kmer_counts = */
  /*   calloc(num_hash_buckets, sizeof(unsigned long)); */
  /* assert(hashed_kmer_counts != NULL); */

  double* pow2 = malloc(num_hyperplanes * sizeof(double));
  assert(pow2 != NULL);
  for (int i = 0; i < num_hyperplanes; ++i) {
    pow2[i] = pow(2, i);
  }

  seq = kseq_init(fp);
  hyperplanes = eigg_hyperplane_ary_init(num_hyperplanes, kmer_len);
  encode_nt   = eigg_encode_ary_init();
  tommy_array_init(hash_bucket_names);
  tommy_hashlin_init(hash_bucket_counts);

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
    if ((num_seqs % 1000) == 0 && num_seqs != 0) {
      fprintf(stderr, "LOG -- reading seqs: %lu\r", num_seqs);
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

    struct Kmers* for_kmers =
      eigg_kmers_new(s1_seq, s1_seq_len, kmer_len);
    double* hashed_for_kmers =
      malloc(for_kmers->num_kmers * sizeof(double));

    /* the rev read might be different lenght */
    struct Kmers* rev_kmers =
      eigg_kmers_new(seq->seq.s, seq->seq.l, kmer_len);
    double* hashed_rev_kmers =
      malloc(rev_kmers->num_kmers * sizeof(double));

    /* hash kmers for forward read */
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

      /* have we seen the hashed_kmer in this sample yet? */
      hash_bucket_count = tommy_hashlin_search(hash_bucket_counts,
                                               hash_bucket_count_compare,
                                               &hashed_kmer,
                                               tommy_inthash_u64(hashed_kmer));

      if (hash_bucket_count) { /* found, increment the count */
        /* TODO worry about overflows? */
        ++hash_bucket_count->count;
      } else { /* not found, add it with a count of 1 */
        ++num_seen_hash_buckets;

        /* track the order */
        /* TODO this could be shared with the name obj in the hash instead? */
        unsigned long* hash_bucket = malloc(sizeof(unsigned long));
        assert(hash_bucket != NULL);
        hash_bucket[0] = hashed_kmer;
        tommy_array_insert(hash_bucket_names, hash_bucket);

        struct hash_bucket_count_t* hb_count =
          malloc(sizeof(struct hash_bucket_count_t));
        assert(hb_count != NULL);

        hb_count->hash_bucket = hashed_kmer;
        hb_count->count = 1;

        tommy_hashlin_insert(hash_bucket_counts,
                             &hb_count->node,
                             hb_count,
                             tommy_inthash_u64(hb_count->hash_bucket));
      }

      /* ++hashed_kmer_counts[(unsigned long)hashed_kmer]; */

      hashed_for_kmers[i] = hashed_kmer;
    }

    /* hash kmers for reverse read */
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

      /* have we seen the hashed_kmer in this sample yet? */
      hash_bucket_count = tommy_hashlin_search(hash_bucket_counts,
                                               hash_bucket_count_compare,
                                               &hashed_kmer,
                                               tommy_inthash_u64(hashed_kmer));

      if (hash_bucket_count) { /* found, increment the count */
        /* TODO worry about overflows? */
        ++hash_bucket_count->count;
      } else { /* not found, add it with a count of 1 */
        ++num_seen_hash_buckets;

        /* track the order */
        /* TODO this could be shared with the name obj in the hash instead? */
        unsigned long* hash_bucket = malloc(sizeof(unsigned long));
        assert(hash_bucket != NULL);
        hash_bucket[0] = hashed_kmer;
        tommy_array_insert(hash_bucket_names, hash_bucket);

        struct hash_bucket_count_t* hb_count =
          malloc(sizeof(struct hash_bucket_count_t));
        assert(hb_count != NULL);

        hb_count->hash_bucket = hashed_kmer;
        hb_count->count = 1;

        tommy_hashlin_insert(hash_bucket_counts,
                             &hb_count->node,
                             hb_count,
                             tommy_inthash_u64(hb_count->hash_bucket));
      }

      /* ++hashed_kmer_counts[(unsigned long)hashed_kmer]; */

      hashed_rev_kmers[i] = hashed_kmer;
    }


    /* each read pair will have a all the kmers of its partner catted
       together with its own */

    /* eigg_print_hashed_record(s1_name, */
    /*                          s1_seq, */
    /*                          s1_qual, */
    /*                          kmer_len, */
    /*                          for_kmers->num_kmers, */
    /*                          rev_kmers->num_kmers, */
    /*                          hashed_for_kmers, */
    /*                          hashed_rev_kmers); */

    /* eigg_print_hashed_record(seq->name.s, */
    /*                          seq->seq.s, */
    /*                          seq->qual.s, */
    /*                          kmer_len, */
    /*                          for_kmers->num_kmers, */
    /*                          rev_kmers->num_kmers, */
    /*                          hashed_for_kmers, */
    /*                          hashed_rev_kmers); */

    /* for (unsigned long i = 0; i < for_kmers->num_kmers; ++i) { */
    /*   fprintf(stdout, */
    /*           "%lu\n", */
    /*           (unsigned long)hashed_for_kmers[i]); */
    /* } */
    /* for (unsigned long i = 0; i < rev_kmers->num_kmers; ++i) { */
    /*   fprintf(stdout, */
    /*           "%lu\n", */
    /*           (unsigned long)hashed_rev_kmers[i]); */
    /* } */

    free(hashed_rev_kmers);
    free(hashed_for_kmers);
    free(s1_name);
    free(s1_seq);
    free(s1_qual);
    eigg_kmers_destroy(rev_kmers);
    eigg_kmers_destroy(for_kmers);
  }
  /* fprintf(stderr, "\n"); */

  /* print kmer counts: outfile is hash_bucket count */
  /* unsigned long count = 0; */

  /* rows, colums, entries */
  fprintf(stdout,
          "1 %lu %lu\n",
          num_hash_buckets,
          num_seen_hash_buckets);
  /* for (unsigned long i = 0; i < num_hash_buckets; ++i) { */
  /*   if ((i % 100000000) == 0 && i != 0 ) { */
  /*     count += 100; */
  /*     fprintf(stderr, */
  /*             "LOG -- processing hash bin counts: %lu million\r", */
  /*             count); */
  /*   } */

  /*   if (hashed_kmer_counts[i] != 0) { */
  /*     fprintf(stdout, */
  /*             "%lu %lu\n", */
  /*             i, */
  /*             hashed_kmer_counts[i]); */
  /*   } */
  /* } */
  /* fprintf(stderr, "\n"); */

  unsigned long* hb_name;
  for (i = 0; i < num_seen_hash_buckets; ++i) {
    if ((i % 100000) == 0 && i != 0) {
      fprintf(stderr,
              "LOG -- processing hash bin counts: %lu\r", i);
    }

    hb_name = (unsigned long*)tommy_array_get(hash_bucket_names, i);
    hash_bucket_count = tommy_hashlin_search(hash_bucket_counts,
                                             hash_bucket_count_compare,
                                             hb_name,
                                             tommy_inthash_u64(hb_name[0]));


    if (hb_name) { /* the hash_bucket is in the hash */
      assert(hb_name[0] == hash_bucket_count->hash_bucket);

      /* sort of like the matrix market format, except the row is not
         printed because we wont know the row until the next step. */
      fprintf(stdout,
              "%lu %lu\n",
              hb_name[0],
              hash_bucket_count->count);
    } else { /* it is missing, but it should be in there */
      fprintf(stderr,
              "Something went wrong! Each bin should be present in bin_counts\n");
      fprintf(stderr,
              "Bin name: %lu\n"
              "Possible hash buckets: %lu\n"
              "Seen hash buckets: %lu\n"
              "Hash bucket num: %lu\n",
              hb_name[0],
              num_hash_buckets,
              num_seen_hash_buckets,
              i);

      return 5;
    }
  }


  kseq_destroy(seq);
  gzclose(fp);
  eigg_hyperplane_ary_destroy(hyperplanes, num_hyperplanes);
  free(encoded_kmer);
  eigg_encode_ary_destroy(encode_nt);
  free(pow2);
  /* free(hashed_kmer_counts); */

  for (tommy_size_t s = 0; s < tommy_array_size(hash_bucket_names); ++s) {
    free(tommy_array_get(hash_bucket_names, s));
  }
  tommy_array_done(hash_bucket_names);
  free(hash_bucket_names);

  tommy_hashlin_foreach(hash_bucket_counts,
                        free);
  tommy_hashlin_done(hash_bucket_counts);
  free(hash_bucket_counts);

  return 0;
}
