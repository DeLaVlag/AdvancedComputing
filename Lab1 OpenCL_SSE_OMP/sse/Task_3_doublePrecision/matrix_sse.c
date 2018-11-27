#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <time.h>
#include <omp.h>

/*****************************************************
the following function generates a "size"-element vector
and a "size x size" matrix
 ****************************************************/
void matrix_vector_gen(int paddedSize, double *matrix, double *vector){

  int i, j=0;

  for(i=0; i<paddedSize; i++)
    vector[i] = ((double)rand())/65535.0; //i*1.2f + 1;//
  for(i=0; i<paddedSize; i++){
    for(j=0; j<paddedSize; j++){
    //   matrix[j] = ((double)rand())/5307.0;
      matrix[i*paddedSize+j] = ((double)rand())/5307.0; //index*1.3f + 1;//      
    }
  }
}

/****************************************************
the following function calculate the below equation
   vector_out = vector_in x matrix_in
 ***************************************************/
void matrix_mult_sq(int size, double *vector_in,
		       double *matrix_in, double *vector_out){
  int rows, cols;
  int j;
  for(cols=0; cols<size; cols++){
    vector_out[cols] = 0.0;
    for(j=0,rows=0; rows<size; j++,rows++)
      vector_out[cols] += vector_in[j] * matrix_in[rows*size + cols];
  }
}

/****************************************************
functie neemt steeds 4 kolommen en loopt 1 voor 1 door de rijen
vermendigvuldigt ze met de volledige vector
telt ze bij elkaar op en stored ze
gaat verder naar de volgende 4 kolommen
 ***************************************************/
void matrix_mult_sse(int size, double *vector_in,
		      double const *matrix_in, double *vector_out){
  __m128d a_line, b_line, r_line;
  int i, j;
  // size = 2 is for double floating point, size = 4 is for single
  for (i=0; i<size; i+=2){
    j = 0;
    // Load 4 values of matrix
    b_line = _mm_load_pd(&matrix_in[i]); // b_line = vec4(matrix[i][0])
    // Broadcast a value to all values in the vector [a,a,a,a];
    a_line = _mm_set1_pd(vector_in[j]);  // a_line = vec4(vector_in[0])
    r_line = _mm_mul_pd(a_line, b_line); // r_line = a_line * b_line
    for (j=1; j<size; j++) {
      b_line = _mm_load_pd(&matrix_in[j*size+i]); // a_line = vec4(column(a, j))
      a_line = _mm_set1_pd(vector_in[j]);  // b_line = vec4(b[i][j])
                                     // r_line += a_line * b_line
      r_line = _mm_add_pd(_mm_mul_pd(a_line, b_line), r_line);
    }
    _mm_store_pd(&vector_out[i], r_line);     // r[i] = r_line
  }
}

// printing the matrix
void print_Matrix(double *matrix, int size){

  int i=0,j;

  printf("\n");
  for (j=0;j<size*size;j++){
    printf(" %1.3f ",matrix[j]);
    i++;
    if (i==size){
      printf("\n");
      i=0;
    }
  }
}

void print_Vector(double *vector, int size){
  int j;
  for (j=0;j<size;j++){
    printf(" %1.3f ",vector[j]);
  }
  printf("\n");
}

int main(int argc, char *argv[]){
  if(argc < 2){
    printf("Usage: %s + matrix size plz\n", argv[0]);
    return 0;
  }

  // int size = atoi(argv[1]);
  // if(size%4 != 0){
  //   printf("This version implements for ""size = 4*n"" only\n");
  //   return 0;
  // }

  int size = atoi(argv[1]);

  int paddedSize = size;

  if (!(size%4 == 0)){
    size = size + (4 - size%4);
  }

  double *vector = (double *)memalign(sizeof(double)*2, sizeof(double)*size);//(double *)malloc(sizeof(double)*size);
  if(vector==NULL){
    printf("can't allocate the required memory for vector\n");
    return 0;
  }

  double *matrix = (double *)memalign(sizeof(double)*2, sizeof(double)*size*size);
  if(matrix==NULL){
    printf("can't allocate the required memory for matrix\n");
    free(vector);
    return 0;
  }

  double *result_sq = (double *)memalign(sizeof(double)*2, sizeof(double)*size);
  if(result_sq==NULL){
    printf("can't allocate the required memory for result_sq\n");
    free(vector);
    free(matrix);
    return 0;
  }

  double *result_pl = (double *)memalign(sizeof(double)*2, sizeof(double)*size);
  if(result_pl==NULL){
    printf("can't allocate the required memory for result_pl\n");
    free(vector);
    free(matrix);
    free(result_sq);
    return 0;
  }

  // init matrix to zero. defaul init gives unknown numbers
  int j;
  for (j=0;j<size*size;j++){
    matrix[j]=0;
  }

  matrix_vector_gen(paddedSize, matrix, vector);

  // print_Matrix(matrix, size);

  double time_sq;
  double time_sse;

  time_sq = omp_get_wtime();
  matrix_mult_sq(size, vector, matrix, result_sq);
  time_sq = omp_get_wtime() - time_sq;

  time_sse = omp_get_wtime();
  matrix_mult_sse(size, vector, matrix, result_pl);
  time_sse = omp_get_wtime() - time_sse;

  // print_Matrix(matrix, size);

  printf("%f\t", time_sq);
  printf("%f\t\n", time_sse);

  //check
  int i;
  for(i=0; i<size; i++)
    if((int)result_sq[i] != (int)result_pl[i]){
      printf("wrong at position %d\n", i);
      // free(vector);
      // free(matrix);
      // free(result_sq);
      // free(result_pl);
      return 0;
    }

  free(vector);
  free(matrix);
  free(result_sq);
  free(result_pl);
  return 1;
}
