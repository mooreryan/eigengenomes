#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "eigg_version.h"

void ary_2d_set_at(unsigned long* ary,
                   unsigned long ridx,
                   unsigned long cidx,
                   unsigned long ncols,
                   unsigned long val)
{
  ary[(ridx * ncols) + cidx] = val;
}

unsigned long ary_2d_get_at(unsigned long* ary,
                            unsigned long ridx,
                            unsigned long cidx,
                            unsigned long ncols)
{
  return ary[(ridx * ncols) + cidx];
}

double l2_norm(unsigned long* ary, unsigned long len)
{
  unsigned long i = 0;
  double norm = 0.0;

  for (i = 0; i < len; ++i) {
    norm += (ary[i] * ary[i]);
  }

  return sqrt(norm);
}

int main(int argc, char *argv[])
{
  if (argc <= 3) {
    fprintf(stderr,
            "VERSION: %s\nUsage: %s "
            "<1: num hyperplanes used in previous step> "
            "*.hash_bucket_counts\n",
            VERSION,
            argv[0]);

    return 1;
  }

  unsigned long hash_bucket = 0;
  unsigned long count = 0;
  unsigned long hash_bucket_count = 0;
  unsigned long i = 0;
  unsigned long j = 0;

  unsigned long val = 0;
  double weighted_count = 0.0;
  double weight = 0.0;
  unsigned long non_zero_samples = 0;
  double hash_bucket_sample_weight = 0.0;


  unsigned int file_idx_offset = 2;
  unsigned long num_samples = argc - file_idx_offset;
  unsigned long num_hash_buckets = pow(2, strtol(argv[1], NULL, 10));

  FILE** infiles = malloc(num_samples * sizeof(FILE*));
  assert(infiles != NULL);

  /* TODO these use way too much memory, consider using a hash table
     for counting instead */
  unsigned long* per_sample_hash_bucket_counts =
    malloc(num_samples * num_hash_buckets * sizeof(unsigned long));

  unsigned long* samples_per_hash_bucket =
    malloc(num_hash_buckets * num_samples * sizeof(unsigned long));

  for (i = 0; i < num_samples * num_hash_buckets; ++i) {
    per_sample_hash_bucket_counts[i] = 0;
    samples_per_hash_bucket[i] = 0;
  }

  double* l2_norms =
    malloc(num_samples * sizeof(double));

  for (i = 0; i < num_samples; ++i) {
    l2_norms[i] = 0.0;
  }

  double* hash_bucket_sample_weights =
    malloc(num_hash_buckets * sizeof(double));

  for (i = 0; i < num_hash_buckets; ++i) {
    hash_bucket_sample_weights[i] = 0.0;
  }

  for (i = file_idx_offset; i < argc; ++i) {
    infiles[i-file_idx_offset] = fopen(argv[i], "r");
    if (!infiles[i-file_idx_offset]) {
      /* TODO handle error */
    }
  }

  for (i = 0; i < num_samples; ++i) {
    fscanf(infiles[i],
           "%lu",
           &num_hash_buckets);

    l2_norms[i] = 0.0;

    /* TODO ensure that num_hash_buckets is the same for all files */

    while (fscanf(infiles[i], "%lu %lu", &hash_bucket, &count) == 2) {
      val = ary_2d_get_at(per_sample_hash_bucket_counts,
                          i,
                          hash_bucket,
                          num_hash_buckets);

      ary_2d_set_at(per_sample_hash_bucket_counts,
                    i,
                    hash_bucket,
                    num_hash_buckets,
                    val + count);

      l2_norms[i] += (count * count);
      ary_2d_set_at(samples_per_hash_bucket,
                    hash_bucket,
                    i,
                    num_samples,
                    1);
    }

    /* see kmer abundance matrix section of the paper */
    l2_norms[i] = sqrt(l2_norms[i]) / sqrt(num_hash_buckets);
  }

  /* see kmer abundance matrix section of the paper */
  for (i = 0; i < num_hash_buckets; ++i) {
    non_zero_samples = 0;
    hash_bucket_sample_weight = 0.0;
    for (j = 0; j < num_samples; ++j) {
      val = ary_2d_get_at(samples_per_hash_bucket, i, j, num_samples);

      if (val > 0) {
        ++non_zero_samples;
      }
    }

    if (non_zero_samples > 0) {

      /* give it a little weight +1 for when the kmers are in all
         samples. TODO this is not in the original...does it
         matter? */
      hash_bucket_sample_weight =
        log2(1 + (num_samples  / (double)non_zero_samples));

    }

    hash_bucket_sample_weights[i] = hash_bucket_sample_weight;
  }

  fprintf(stdout,
          "%lu %lu\n",
          num_samples,
          num_hash_buckets);
  for (i = 0; i < num_samples; ++i) {
    for (j = 0; j < num_hash_buckets; ++j) {
      /* log normalization of tf count... TOOD not in original, better
         of worse? */
      /* hash_bucket_count = */
      /*   1 + log2(ary_2d_get_at(per_sample_hash_bucket_counts, */
      /*                          i, */
      /*                          j, */
      /*                          num_hash_buckets)); */

      hash_bucket_count = ary_2d_get_at(per_sample_hash_bucket_counts,
                                        i,
                                        j,
                                        num_hash_buckets);

      if (hash_bucket_count != 0) {
        weight = hash_bucket_sample_weights[j] / l2_norms[i];
        weighted_count = hash_bucket_count * weight;
        fprintf(stdout,
                "%lu %lu %.5f\n",
                i,
                j,
                weighted_count);
      }
    }
  }

  for (i = 0; i < num_samples; ++i) {
    fclose(infiles[i]);
  }
  free(infiles);
  free(per_sample_hash_bucket_counts);
  free(samples_per_hash_bucket);
  free(hash_bucket_sample_weights);
  free(l2_norms);

  return 0;
}
