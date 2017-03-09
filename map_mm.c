#include <stdio.h>

int main(int argc, char *argv[])
{
  unsigned long actual_num_hbs = 0;
  unsigned long i = 0;
  unsigned long hb_name = 0;
  unsigned long num_entries = 0;
  unsigned long num_hbs = 0;
  unsigned long num_samples = 0;
  unsigned long sample_name = 0;

  /* hb_names start at 1 */
  unsigned long prev_hb_name = 0;
  unsigned long mapped_hb_name = 0;


  double count = 0.0;

  fscanf(stdin,
         "%lu %lu %lu %lu",
         &actual_num_hbs,
         &num_hbs,
         &num_samples,
         &num_entries);

  fprintf(stderr,
          "%lu\n",
          actual_num_hbs);

  fprintf(stdout,
          "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf(stdout,
          "%lu %lu %lu\n",
          actual_num_hbs,
          num_samples,
          num_entries);


  for (i = 0; i < num_entries; ++i) {
    fscanf(stdin,
           "%lu %lu %lf",
           &hb_name,
           &sample_name,
           &count);

    if (hb_name == prev_hb_name) {
      fprintf(stdout,
              "%lu %lu %.5lf\n",
              mapped_hb_name,
              sample_name,
              count);
    } else {

      ++mapped_hb_name;

      fprintf(stderr,
              "%lu %lu\n",
              mapped_hb_name,
              hb_name);

      fprintf(stdout,
              "%lu %lu %.5lf\n",
              mapped_hb_name,
              sample_name,
              count);
    }
    prev_hb_name = hb_name;
  }

  return 0;
}
