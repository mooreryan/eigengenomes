#ifndef EIGG_HYPERPLANE_H
#define EIGG_HYPERPLANE_H

double** eigg_hyperplane_ary_init(int num_hyperplanes, int num_coef);
void eigg_hyperplane_ary_destroy(double** hyperplane_ary, int num_hyperplanes);

#endif
