#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "eigg_version.h"
#include "vendor/tommyarray.h"
#include "vendor/tommyhashlin.h"

struct item_t {
  unsigned long key;
  double val;

  tommy_node node;
};

int item_compare(const void* arg, const void* item)
{
  return *(const unsigned long*)arg !=
    ((const struct item_t*)item)->key;
}

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

  /* to hold the tommy_hashlin_search() return vals */
  struct item_t* item = malloc(sizeof(struct item_t));

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
  unsigned long num_hash_buckets =
    (unsigned long)round(pow(2, strtol(argv[1], NULL, 10)));

  FILE** infiles = malloc(num_samples * sizeof(FILE*));
  assert(infiles != NULL);

  /* TODO these use way too much memory, consider using a hash table
     for counting instead */
  /* unsigned long* per_sample_hash_bucket_counts = */
  /*   malloc(num_samples * num_hash_buckets * sizeof(unsigned long)); */

  tommy_hashlin** per_sample_hash_bucket_counts
    = malloc(num_samples * sizeof(tommy_hashlin*));
  assert(per_sample_hash_bucket_counts != NULL);
  for (i = 0; i < num_samples; ++i) {
    per_sample_hash_bucket_counts[i]
      = malloc(sizeof(tommy_hashlin));

    tommy_hashlin_init(per_sample_hash_bucket_counts[i]);
  }

  unsigned long* samples_per_hash_bucket =
    malloc(num_hash_buckets * num_samples * sizeof(unsigned long));

  for (i = 0; i < num_samples * num_hash_buckets; ++i) {
    /* per_sample_hash_bucket_counts[i] = 0; */
    samples_per_hash_bucket[i] = 0;
  }

  double* l2_norms =
    malloc(num_samples * sizeof(double));

  for (i = 0; i < num_samples; ++i) {
    l2_norms[i] = 0.0;
  }


  tommy_array* hash_bucket_names
    = malloc(sizeof(tommy_array));
  assert(hash_bucket_names != NULL);
  tommy_array_init(hash_bucket_names);

  tommy_hashlin* hash_bucket_sample_weights
    = malloc(sizeof(tommy_hashlin));
  assert(hash_bucket_sample_weights != NULL);
  tommy_hashlin_init(hash_bucket_sample_weights);

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

    /* i is sample idx */
    while (fscanf(infiles[i], "%lu %lu", &hash_bucket, &count) == 2) {
      /* TODO check and see if this sample's hash bucket name has been
         seen before */
      struct item_t* new_item = malloc(sizeof(struct item_t));
      new_item->key = hash_bucket;
      new_item->val = count;
      tommy_hashlin_insert(per_sample_hash_bucket_counts[i],
                           &new_item->node,
                           new_item,
                           tommy_inthash_u64(new_item->key));

      /* val = ary_2d_get_at(per_sample_hash_bucket_counts, */
      /*                     i, */
      /*                     hash_bucket, */
      /*                     num_hash_buckets); */

      /* ary_2d_set_at(per_sample_hash_bucket_counts, */
      /*               i, */
      /*               hash_bucket, */
      /*               num_hash_buckets, */
      /*               val + count); */

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

      /* insert it into the tracking array */
      unsigned long* hash_bucket_name = malloc(sizeof(unsigned long));
      assert(hash_bucket_name != NULL);
      hash_bucket_name[0] = i;
      tommy_array_insert(hash_bucket_names, hash_bucket_name);

      /* insert it into the hash table */
      struct item_t* new_item = malloc(sizeof(struct item_t));
      new_item->key = i;
      new_item->val = hash_bucket_sample_weight;

      tommy_hashlin_insert(hash_bucket_sample_weights,
                           &new_item->node,
                           new_item,
                           tommy_inthash_u64(new_item->key));
    }

    /* hash_bucket_sample_weights[i] = hash_bucket_sample_weight; */
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

      item = tommy_hashlin_search(per_sample_hash_bucket_counts[i],
                                  item_compare,
                                  &j,
                                  tommy_inthash_u64(j));
      if (!item) {
        fprintf(stderr,
                "WARN -- hash bucket %lu was not found in pshbc hash\n",
                j);

        /* return 5; */
      } else {
        hash_bucket_count = (unsigned long)round(item->val);

        /* hash_bucket_count = ary_2d_get_at(per_sample_hash_bucket_counts, */
        /*                                   i, */
        /*                                   j, */
        /*                                   num_hash_buckets); */

        if (hash_bucket_count != 0) {
          /* get the sample weight */
          item = tommy_hashlin_search(hash_bucket_sample_weights,
                                      item_compare,
                                      &j,
                                      tommy_inthash_u64(j));
          if (!item) {
            fprintf(stderr,
                    "ERROR -- hash bucket %lu was not found in weight hash\n",
                    j);

            return 5;
          } else {
            /* TODO there is something wonky in valgrind related to
               this */
            hash_bucket_sample_weight = item->val;
          }

          weight = hash_bucket_sample_weight / l2_norms[i];
          weighted_count = hash_bucket_count * weight;
          fprintf(stdout,
                  "%lu %lu %.5f\n",
                  i,
                  j,
                  weighted_count);
        }
      }
    }
  }

  for (i = 0; i < num_samples; ++i) {
    fclose(infiles[i]);
  }
  free(infiles);
  free(samples_per_hash_bucket);
  free(l2_norms);

  tommy_array_done(hash_bucket_names);
  free(hash_bucket_names);
  tommy_hashlin_done(hash_bucket_sample_weights);
  free(hash_bucket_sample_weights);
  free(item);

  for (i = 0; i < num_samples; ++i) {
    tommy_hashlin_done(per_sample_hash_bucket_counts[i]);
    free(per_sample_hash_bucket_counts[i]);

    /* TODO free the elements of the hash? */
  }
  free(per_sample_hash_bucket_counts);

  /* TODO free the items in hash_bucket_names and hash_bucket_sample_weights */

  return 0;
}
