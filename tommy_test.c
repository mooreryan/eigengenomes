#include <stdlib.h>
#include <stdio.h>

#include "vendor/tommyarray.h"
#include "vendor/tommyhashlin.h"

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

int main(int argc, char *argv[])
{

  tommy_array*   hash_bucket_names = malloc(sizeof(tommy_array));
  assert(hash_bucket_names != NULL);

  tommy_hashlin* hash_bucket_counts = malloc(sizeof(tommy_hashlin));
  assert(hash_bucket_counts != NULL);

  unsigned long i = 0;

  struct hash_bucket_count_t* hb_count =
    malloc(sizeof(struct hash_bucket_count_t));
  assert(hb_count != NULL);

  tommy_array_init(hash_bucket_names);
  tommy_hashlin_init(hash_bucket_counts);

  hb_count = tommy_hashlin_search(hash_bucket_counts,
                                  hash_bucket_count_compare,
                                  &i,
                                  tommy_inthash_u64(i));

  if (hb_count) {
    printf("I found %lu\n", i);
  } else {
    printf("I did NOT find %lu\n", i);
    printf("I will insert it!\n");

    struct hash_bucket_count_t* new_hb_count =
      malloc(sizeof(struct hash_bucket_count_t));
    assert(new_hb_count != NULL);

    new_hb_count->hash_bucket = i;
    new_hb_count->count = 1;

    tommy_hashlin_insert(hash_bucket_counts,
                         &new_hb_count->node,
                         new_hb_count,
                         tommy_inthash_u64(new_hb_count->hash_bucket));
  }

  hb_count = tommy_hashlin_search(hash_bucket_counts,
                                  hash_bucket_count_compare,
                                  &i,
                                  tommy_inthash_u64(i));

  puts("");
  if (hb_count) {
    printf("I found %lu, with count: %lu\n",
           hb_count->hash_bucket,
           hb_count->count);
    printf("incrementing the count\n");
    ++hb_count->count;
    printf("count is now: %lu\n",
           hb_count->count);
  } else {
    printf("I did NOT find %lu\n", i);
  }

  hb_count = tommy_hashlin_search(hash_bucket_counts,
                                  hash_bucket_count_compare,
                                  &i,
                                  tommy_inthash_u64(i));

  puts("");
  if (hb_count) {
    printf("I found %lu, with count: %lu\n",
           hb_count->hash_bucket,
           hb_count->count);
    printf("incrementing the count\n");
    ++hb_count->count;
    printf("count is now: %lu\n",
           hb_count->count);
  } else {
    printf("I did NOT find %lu\n", i);
  }


  tommy_array_done(hash_bucket_names);
  free(hash_bucket_names);
  tommy_hashlin_done(hash_bucket_counts);
  free(hash_bucket_counts);
  free(hb_count);
  return 0;
}
