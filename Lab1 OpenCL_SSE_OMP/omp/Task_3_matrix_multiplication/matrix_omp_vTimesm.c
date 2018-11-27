/*
This code is compiled and executed via a couple of shell scripts
resided in the folder of this file.

The code in this file is almost the same as the original code.
Changed code can be recognized by the comments which start with 3
backslashes (///).

The code has to be executed via a terminal and typing: sh ExecuteMe.sh
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/*****************************************************
the following function generates a "size"-element vector
and a "size x size" matrix
 ****************************************************/
void matrix_vector_gen(int size, double *matrix, double *vector) {
  int i;
  for (i = 0; i < size; i++)
    vector[i] = ((double)rand()) / 65535.0;
  for (i = 0; i < size * size; i++)
    matrix[i] = ((double)rand()) / 5307.0;
}

/****************************************************
the following function calculates the below equation
     vector_out = vector_in x matrix_in
 ***************************************************/
void matrix_mult_sq(int size, double *vector_in,
                    double *matrix_in, double *vector_out) {
  int rows, cols;
  int j;
  for (cols = 0; cols < size; cols++) {
    vector_out[cols] = 0.0;
    for (j = 0, rows = 0; rows < size; j++, rows++)
      vector_out[cols] += vector_in[j] * matrix_in[rows * size + cols];
  }
}

/****************************************************
the following function calculates the equation below
     vector_out = vector_in x matrix_in
     With the OMP parallellization.
 ***************************************************/
void matrix_mult_pl(int size, double *vector_in,
                    double *matrix_in, double *vector_out) {
  int rows, cols;
  int j;
  /* shared() the input parameters and output variabel*/
  /* privat() the local variables */
  # pragma omp parallel             \
  shared(size, vector_in, matrix_in, vector_out)    \
  private(rows, cols, j)
  # pragma omp for
  for (cols = 0; cols < size; cols++) {
    vector_out[cols] = 0.0;
    for (j = 0, rows = 0; rows < size; j++, rows++)
      vector_out[cols] += vector_in[j] * matrix_in[rows * size + cols];
  }
}

/****************************************************
the main function.
 ***************************************************/
int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("Usage: %s matrix/vector_size\n", argv[0]);
    return 0;
  }

  int size = atoi(argv[1]);
  double *vector    = (double *)malloc(sizeof(double) * size);
  double *matrix    = (double *)malloc(sizeof(double) * size * size);
  double *result_sq = (double *)malloc(sizeof(double) * size);
  double *result_pl = (double *)malloc(sizeof(double) * size);
  matrix_vector_gen(size, matrix, vector);

  double time_sq = 0;
  double time_pl = 0;

  time_sq = omp_get_wtime();
  matrix_mult_sq(size, vector, matrix, result_sq);
  time_sq = omp_get_wtime() - time_sq;

  time_pl = omp_get_wtime();
  matrix_mult_pl(size, vector, matrix, result_pl);
  time_pl = omp_get_wtime() - time_pl;

  ///Calculate the performance differences between normal execution and 
  ///This section does not influence the performance measurement.
  printf("SEQ: \t    %f (sec)\t", time_sq);
  printf("PAR: \t(t:%d) \t (p:%d) \t %f (sec)\t",
         omp_get_max_threads(), omp_get_num_procs(), time_pl);
  printf("SEQ / PAR: \t    %f ", time_sq / time_pl);

  //check, the calculation SER. vs. PAR.
  int i;
  for (i = 0; i < size; i++) {
    if (result_sq[i] != result_pl[i]) {
      printf("wrong at position %d\n", i);
      return 0;
    }
  }

  free(vector);
  free(matrix);
  free(result_sq);
  free(result_pl);
  return 1;
}
