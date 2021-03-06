#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "eigg_version.h"
#include "vendor/tommyhashlin.h"
#include "vendor/tommylist.h"

struct sphb_item_t {
  unsigned long hash_bucket;
  unsigned int* samples;
  unsigned long num_samples;

  tommy_node node;
};

struct sphb_item_t*
sphb_item_create(unsigned long hash_bucket,
                 unsigned long num_samples)
{
  struct sphb_item_t* sphb_item =
    malloc(sizeof(struct sphb_item_t));
  assert(sphb_item != NULL);

  sphb_item->hash_bucket = hash_bucket;

  unsigned int* samples =
    malloc(num_samples * sizeof(int));

  for (unsigned long i = 0; i < num_samples; ++i) {
    samples[i] = 0;
  }

  sphb_item->samples = samples;
  sphb_item->num_samples = num_samples;

  return sphb_item;
}

void
sphb_item_destroy(struct sphb_item_t* sphb_item)
{
  free(sphb_item->samples);
  free(sphb_item);
}

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

struct lu_t {
  long unsigned hash_bucket_name;
  tommy_node node;
};

int lu_sort_func(const void* a, const void* b)
{
  unsigned long bucket_a = 0;
  unsigned long bucket_b = 0;

  bucket_a = ((struct lu_t*)a)->hash_bucket_name;
  bucket_b = ((struct lu_t*)b)->hash_bucket_name;

  if (bucket_a < bucket_b) {
    return -1;
  } else if (bucket_a == bucket_b) {
    return 0;
  } else {
    return 1;
  }
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
  if (argc <= 2) {
    fprintf(stderr,
            "VERSION: %s\nUsage: %s "
            "<1: num hyperplanes used in previous step> "
            "*.hash_bucket_counts\n",
            VERSION,
            argv[0]);

    return 1;
  }

  /* to hold the tommy_hashlin_search() return vals */
  struct item_t* item;

  struct sphb_item_t* sphb_item;

  unsigned long rows = 0;
  unsigned long cols = 0;
  unsigned long entries = 0;
  unsigned long total_entries = 0;


  unsigned long hash_bucket = 0;
  unsigned long count = 0;
  unsigned long hash_bucket_count = 0;
  unsigned long i = 0;
  unsigned long line = 0;
  unsigned long s_i = 0;
  unsigned long hb_i = 0;

  struct lu_t* hash_bucket_name;

  unsigned long actual_num_hash_buckets = 0;
  unsigned long s = 0;

  tommy_node* thingy;

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

  /* each sample hash its own counting hash to track the number of
     kmers assigned to each hash bucket */
  tommy_hashlin** per_sample_hash_bucket_counts
    = malloc(num_samples * sizeof(tommy_hashlin*));
  assert(per_sample_hash_bucket_counts != NULL);
  for (s_i = 0; s_i < num_samples; ++s_i) {
    per_sample_hash_bucket_counts[s_i]
      = malloc(sizeof(tommy_hashlin));
    assert(per_sample_hash_bucket_counts[s_i] != NULL);

    tommy_hashlin_init(per_sample_hash_bucket_counts[s_i]);
  }

  /* for each hash bucket, track how many samples had a least one kmer
     in this hash bucket */
  tommy_hashlin* samples_per_hash_bucket =
    malloc(sizeof(tommy_hashlin));
  assert(samples_per_hash_bucket != NULL);
  tommy_hashlin_init(samples_per_hash_bucket);

  double* l2_norms =
    malloc(num_samples * sizeof(double));

  for (s_i = 0; s_i < num_samples; ++s_i) {
    l2_norms[s_i] = 0.0;
  }

  tommy_list* hash_bucket_names
    = malloc(sizeof(tommy_list));
  assert(hash_bucket_names != NULL);
  tommy_list_init(hash_bucket_names);

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

  for (s_i = 0; s_i < num_samples; ++s_i) {
    fscanf(infiles[s_i],
           "%lu %lu %lu",
           &rows, &cols, &entries);
    total_entries += entries;

    num_hash_buckets = cols;

    l2_norms[s_i] = 0.0;

    /* TODO ensure that num_hash_buckets is the same for all files */

    /* read all lines in this sample file */
    line = 0;
    for (line = 0; line < entries; ++line) {

      /* TODO check that the correct number of items are read */
      fscanf(infiles[s_i], "%lu %lu", &hash_bucket, &count);
      if ((line % 1000000) == 0) {
        fprintf(stderr,
                "Reading and tracking sample %lu of %lu = %.2f%%, line #%lu\r",
                s_i+1,
                num_samples,
                (s_i+1) / (double)num_samples * 100,
                line);
      }

      /* TODO check and see if this sample's hash bucket name has been
         seen before */
      struct item_t* new_item = malloc(sizeof(struct item_t));
      new_item->key = hash_bucket;
      new_item->val = count;
      tommy_hashlin_insert(per_sample_hash_bucket_counts[s_i],
                           &new_item->node,
                           new_item,
                           tommy_inthash_u64(new_item->key));

      /* val = ary_2d_get_at(per_sample_hash_bucket_counts, */
      /*                     s_i, */
      /*                     hash_bucket, */
      /*                     num_hash_buckets); */

      /* ary_2d_set_at(per_sample_hash_bucket_counts, */
      /*               s_i, */
      /*               hash_bucket, */
      /*               num_hash_buckets, */
      /*               val + count); */

      l2_norms[s_i] += (count * count);

      /* TODO could this be combined with the previous hash? */

      sphb_item = tommy_hashlin_search(samples_per_hash_bucket,
                                       item_compare,
                                       &hash_bucket,
                                       tommy_inthash_u64(hash_bucket));
      if (sphb_item) { /* has been seen */
        /* this hashbucket has been seen in this sample */
        sphb_item->samples[s_i] = 1;
      } else {
        /* insert hash bucket into the tracking list */
        hash_bucket_name = malloc(sizeof(struct lu_t));
        assert(hash_bucket_name != NULL);
        hash_bucket_name->hash_bucket_name = hash_bucket;
        ++actual_num_hash_buckets;
        tommy_list_insert_head(hash_bucket_names, &hash_bucket_name->node, hash_bucket_name);


        struct sphb_item_t* new_sphb_item =
          sphb_item_create(hash_bucket, num_samples);

        /* this hashbucket has been seen in this sample */
        new_sphb_item->samples[s_i] = 1;

        tommy_hashlin_insert(samples_per_hash_bucket,
                             &new_sphb_item->node,
                             new_sphb_item,
                             tommy_inthash_u64(new_sphb_item->hash_bucket));
      }
    }

    /* see kmer abundance matrix section of the paper */
    l2_norms[s_i] = sqrt(l2_norms[s_i]) / sqrt(num_hash_buckets);
  }
  fprintf(stderr, "\n");

  /* see kmer abundance matrix section of the paper */
  fprintf(stderr,
          "INFO -- possible hash buckets: %lu, actual hash buckets: %lu, "
          "%.5f%% of possible\n",
          num_hash_buckets,
          actual_num_hash_buckets,
          actual_num_hash_buckets / (double)num_hash_buckets * 100);

  thingy = tommy_list_head(hash_bucket_names);
  s = 0;
  while (thingy) {
    ++s;
    struct lu_t* lu = thingy->data; /* get the obj pointer */
    hb_i = lu->hash_bucket_name;
    thingy = thingy->next; /* go to the next element */

    if ((s % 1000000) == 0) {
      fprintf(stderr,
              "Global weighting counts -- %lu of %lu = %.2f%%\r",
              s,
              actual_num_hash_buckets,
              s / (double) actual_num_hash_buckets * 100);
    }

    non_zero_samples = 0;
    hash_bucket_sample_weight = 0.0;
    for (s_i = 0; s_i < num_samples; ++s_i) {
      /* val = ary_2d_get_at(samples_per_hash_bucket, hb_i, s_i, num_samples); */

      sphb_item = tommy_hashlin_search(samples_per_hash_bucket,
                                       item_compare,
                                       &hb_i,
                                       tommy_inthash_u64(hb_i));

      if (!sphb_item) {
        /* fprintf(stderr, */
        /*         "WARN -- hash bucket %lu was not found in sphb hash\n", */
        /*         hb_i); */
      } else {
        if (sphb_item->samples[s_i] == 1) {
          ++non_zero_samples;
        }
      }
    }

    if (non_zero_samples > 0) {

      /* give it a little weight +1 for when the kmers are in all
         samples. TODO this is not in the original...does it
         matter? */
      hash_bucket_sample_weight =
        log2(1 + (num_samples  / (double)non_zero_samples));

      /* insert it into the hash table */
      struct item_t* new_item = malloc(sizeof(struct item_t));
      new_item->key = hb_i;
      new_item->val = hash_bucket_sample_weight;

      tommy_hashlin_insert(hash_bucket_sample_weights,
                           &new_item->node,
                           new_item,
                           tommy_inthash_u64(new_item->key));
    }

    /* hash_bucket_sample_weights[hb_i] = hash_bucket_sample_weight; */
  }
  fprintf(stderr,
          "Global weighting counts -- %lu of %lu = %.2f%%\n",
          s,
          actual_num_hash_buckets,
          s / (double) actual_num_hash_buckets * 100);

  /* fprintf(stdout, "%%%%MatrixMarket matrix coordinate real general\n"); */
  fprintf(stdout,
          "%lu %lu %lu %lu\n",
          actual_num_hash_buckets,
          num_hash_buckets,
          num_samples,
          total_entries);


  /* sort the output in hash bucket order. TODO could try making a map
     from hash bucket name down to 0..actual_num_hash_buckets */
  fprintf(stderr,
          "Sorting hash bucket names\n");
  tommy_list_sort(hash_bucket_names, lu_sort_func);

  thingy = tommy_list_head(hash_bucket_names);

  /* print out info and do final weighting */
  s = 0;
  while (thingy) {
    ++s;
    struct lu_t* lu = thingy->data; /* get the obj pointer */
    hb_i = lu->hash_bucket_name;
    thingy = thingy->next; /* go to the next element */

    if ((s % 1000000) == 0) {
      fprintf(stderr,
              "Final weighting counts, hash bucket %lu of %lu = %.2f%%\r",
              s,
              actual_num_hash_buckets,
              s / (double) actual_num_hash_buckets * 100);
    }

    for (s_i = 0; s_i < num_samples; ++s_i) {

      /* log normalization of tf count... TOOD not in original, better
         of worse? */
      /* hash_bucket_count = */
      /*   1 + log2(ary_2d_get_at(per_sample_hash_bucket_counts, */
      /*                          s_i, */
      /*                          hb_i, */
      /*                          num_hash_buckets)); */

      item = tommy_hashlin_search(per_sample_hash_bucket_counts[s_i],
                                  item_compare,
                                  &hb_i,
                                  tommy_inthash_u64(hb_i));
      if (!item) {
        /* fprintf(stderr, */
        /*         "WARN -- hash bucket %lu was not found in pshbc hash\n", */
        /*         hb_i); */

        /* return 5; */
      } else {
        hash_bucket_count = (unsigned long)round(item->val);

        /* hash_bucket_count = ary_2d_get_at(per_sample_hash_bucket_counts, */
        /*                                   s_i, */
        /*                                   hb_i, */
        /*                                   num_hash_buckets); */

        if (hash_bucket_count != 0) {
          /* get the sample weight */
          item = tommy_hashlin_search(hash_bucket_sample_weights,
                                      item_compare,
                                      &hb_i,
                                      tommy_inthash_u64(hb_i));
          if (!item) {
            fprintf(stderr,
                    "ERROR -- hash bucket %lu was not found in weight hash\n",
                    hb_i);

            return 5;
          } else {
            /* TODO there is something wonky in valgrind related to
               this */
            hash_bucket_sample_weight = item->val;
          }

          weight = hash_bucket_sample_weight / l2_norms[s_i];
          weighted_count = hash_bucket_count * weight;
          /* the MatrixMarket format has 1-based indices...ie A[1,1]
             is the first entry in Matrix A */
          fprintf(stdout,
                  "%lu %lu %.5f\n",
                  hb_i+1,
                  s_i+1,
                  weighted_count);
        }
      }
    }
  }
  fprintf(stderr,
          "Final weighting counts, hash bucket %lu of %lu = %.2f%%\n",
          s,
          actual_num_hash_buckets,
          s / (double) actual_num_hash_buckets * 100);


  fprintf(stderr, "Cleaning up!\n");

  for (s_i = 0; s_i < num_samples; ++s_i) {
    fclose(infiles[s_i]);
  }
  free(infiles);

  tommy_hashlin_foreach(samples_per_hash_bucket,
                        (tommy_foreach_func*)sphb_item_destroy);
  tommy_hashlin_done(samples_per_hash_bucket);
  free(samples_per_hash_bucket);

  free(l2_norms);

  tommy_hashlin_foreach(hash_bucket_sample_weights,
                        free);
  tommy_hashlin_done(hash_bucket_sample_weights);
  free(hash_bucket_sample_weights);

  for (s_i = 0; s_i < num_samples; ++s_i) {
    tommy_hashlin_foreach(per_sample_hash_bucket_counts[s_i],
                          free); /* APPLE is this the correct free func? */
    tommy_hashlin_done(per_sample_hash_bucket_counts[s_i]);
    free(per_sample_hash_bucket_counts[s_i]);

    /* TODO free the elements of the hash? */
  }
  free(per_sample_hash_bucket_counts);

  tommy_list_foreach(hash_bucket_names, free);
  free(hash_bucket_names);

  fprintf(stderr, "Finished %s!\n", argv[0]);
  return 0;
}
