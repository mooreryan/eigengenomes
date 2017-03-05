#include "eigg_hyperplane.h"

#include <stdlib.h>

static double
rand_real()
{
  /* TODO switch to random(4) for better random-ness? */
  return (double) rand() / (double) RAND_MAX;
}

/* random hyperplane through the origin */
static double*
hyperplane_init(int num_coef)
{
  double* hyperplane = malloc(num_coef * sizeof(double));

  int i = 0;
  double coef = 0.0;

  for (i = 0; i < num_coef; ++i) {
    coef = rand_real();

    /* if (rand_real() > 0.5) { */
    /*   coef = 1 / coef; */
    /* } */

    if (rand_real() > 0.5) {
      coef = -coef;
    }

    hyperplane[i] = coef;
  }

  return hyperplane;
}

static void
hyperplane_destroy(double* hyperplane)
{
  free(hyperplane);
}

/* num_coef is normally the kmer len */
double**
eigg_hyperplane_ary_init(int num_hyperplanes,
                         int num_coef)
{
  int i = 0;

  double** hyperplane_ary = malloc(num_hyperplanes * sizeof(double*));

  for (i = 0; i < num_hyperplanes; ++i) {
    hyperplane_ary[i] = hyperplane_init(num_coef);
  }

  return hyperplane_ary;
}

void
eigg_hyperplane_ary_destroy(double** hyperplane_ary,
                            int num_hyperplanes)
{
  int i = 0;

  for (i = 0; i < num_hyperplanes; ++i) {
    hyperplane_destroy(hyperplane_ary[i]);
  }

  free(hyperplane_ary);
}
