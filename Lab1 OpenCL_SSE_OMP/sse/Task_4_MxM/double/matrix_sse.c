#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <time.h>
#include <omp.h>

/*****************************************************
the following function generates a "size"-element vector
and a "size x size" matrix
 ****************************************************/
void matrix_vector_gen(int paddedSize, double *matrix){

  int i,j;

  for(i=0; i<paddedSize; i++){
    for(j=0; j<paddedSize; j++){
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
function takes 4 values from first 4 columns 
multiply 4 elements with one scalar of 4 elements
add them en stores them
moves to the next 4 colomns and iterates over it
 ***************************************************/
void matrix_mult_sse(int size, double *matrix_Neo, double *matrix_Tri, double *matrix_Mrph){
  __m128d a_line, b_line, r_line;
  int i, j, k, index=0;
  // size = 2 is for double, size = 4 is for float
  for (i=0; i<size*size; i+=2){
    // Load 4 values of matrix
    a_line = _mm_load_pd(matrix_Neo);
    b_line = _mm_set1_pd(matrix_Tri[i]);
    // Broadcast a value to all values in the vector [a,a,a,a];
    r_line = _mm_mul_pd(a_line, b_line);
    for (j=1; j<size; j++) {
      a_line = _mm_load_pd(&matrix_Neo[j*size]);
      b_line = _mm_set1_pd(matrix_Tri[i+j]);
      r_line = _mm_add_pd(_mm_mul_pd(a_line, b_line), r_line);
    }
    // this function stores the results one after another in memory
    _mm_store_pd(&matrix_Mrph[i], r_line);
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

// init matrix to zero. defaul init gives unknown numbers
void init_Matrix(double *matrix, int size){
  int j;
  for (j=0;j<size*size;j++){
    matrix[j]=0;
  }
}

int main(int argc, char *argv[]){
  if(argc < 2){
    printf("Usage: %s + matrix size plz\n", argv[0]);
    return 0;
  }

  int size = atoi(argv[1]);

  int paddedSize = size;

  // Checking size for padding
  if (!(size%4 == 0)){
    size = size + (4 - size%4);
  }

  double *vector = (double *)memalign(sizeof(double)*2, sizeof(double)*size);//(double *)malloc(sizeof(double)*size);
  if(vector==NULL){
    printf("can't allocate the required memory for vector\n");
    return 0;
  }

  double *matrix_Neo = (double *)memalign(sizeof(double)*2, sizeof(double)*size*size);
  if(matrix_Neo==NULL){
    printf("can't allocate the required memory for matrix\n");
    return 0;
  }

  double *matrix_Tri = (double *)memalign(sizeof(double)*2, sizeof(double)*size*size);
  if(matrix_Tri==NULL){
    printf("can't allocate the required memory for matrix\n");
    free(matrix_Neo);
    return 0;
  }

  double *matrix_Mrph = (double *)memalign(sizeof(double)*2, sizeof(double)*size*size*size);
  if(matrix_Mrph==NULL){
    printf("can't allocate the required memory for matrix\n");
    free(matrix_Tri);
    free(matrix_Neo);
    return 0;
  }

  // init_Matrix(matrix, size);
  init_Matrix(matrix_Neo, size);
  init_Matrix(matrix_Tri, size);
  init_Matrix(matrix_Mrph, size);

  // matrix_vector_gen(size, matrix, vector, rows, cols);
  matrix_vector_gen(paddedSize, matrix_Neo);
  matrix_vector_gen(paddedSize, matrix_Tri);

  double time_sq;
  double time_sse;

  // time_sq = omp_get_wtime();
  // matrix_mult_sq(size, vector, matrix, result_sq);
  // time_sq = omp_get_wtime() - time_sq;

  time_sse = omp_get_wtime();
  matrix_mult_sse(size, matrix_Neo, matrix_Tri, matrix_Mrph);
  time_sse = omp_get_wtime() - time_sse;

  print_Matrix(matrix_Neo, size);
  print_Matrix(matrix_Tri, size);
  print_Matrix(matrix_Mrph, size);

  // printf("SEQ: %f (sec)\t", time_sq);
  // printf("PAR: %f (sec)\t", time_sse);
  // printf("SEQ / PAR: %f \n", time_sq/time_sse);

  //check
  /*int i;
  for(i=0; i<size; i++)
    if((int)result_sq[i] != (int)result_pl[i]){
      printf("wrong at position %d\n", i);
      free(vector);
      free(matrix);
      free(result_sq);
      free(result_pl);
      return 0;
    }*/

  free(vector);
  free(matrix_Neo);
  free(matrix_Tri);
  free(matrix_Mrph);
  return 1;
}
