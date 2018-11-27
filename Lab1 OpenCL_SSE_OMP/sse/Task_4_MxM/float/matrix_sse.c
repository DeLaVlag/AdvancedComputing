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
void matrix_vector_gen(int size, float *matrix_in){

  int i,j, index=0;

  for(i=0; i<size*size; i++)
  matrix_in[i] = i*1.3f + 1;//((float)rand())/5307.0f;

}

/****************************************************
the following function calculate the below equation
   vector_out = vector_in x matrix_in
 ***************************************************/
void matrix_mult_sq(int size, float *vector_in,
		       float *matrix_in, float *vector_out){
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
void matrix_mult_sse(int size, float *matrix_1, float *matrix_2, float *matrix_3){
  __m128 a_line, b_line, r_line;
  int i, j, k, index=0;
  for (i=0; i<size*size; i+=4){
    a_line = _mm_load_ps(matrix_1); 
    b_line = _mm_set1_ps(matrix_2[i]);
    r_line = _mm_mul_ps(a_line, b_line);
    for (j=1; j<size; j++) {
      a_line = _mm_load_ps(&matrix_1[j*size]);
      b_line = _mm_set1_ps(matrix_2[i+j]);
      r_line = _mm_add_ps(_mm_mul_ps(a_line, b_line), r_line);
    }
    _mm_store_ps(&matrix_3[i], r_line);
  }
}

// printing the matrix
void print_Matrix(float *matrix_in, int size){

  int i=0,j;

  printf("\n");
  for (j=0;j<size*size;j++){
    printf(" %1.3f ",matrix_in[j]);
    i++;
    if (i==size){
      printf("\n");
      i=0;
    }
  }
}

// init matrix to zero. default init gives unknown numbers
void init_Matrix(float *matrix_in, int size){
  int j;
  for (j=0;j<size*size;j++){
    matrix_in[j]=0;
  }
}

int main(int argc, char *argv[]){
  if(argc < 2){
    printf("Usage: %s Enter matrix size plz\n", argv[0]);
    return 0;
  }

  int size = atoi(argv[1]);
  if(size%4 != 0){
    printf("This version implements for ""size = 4*n"" only\n");
    return 0;
  }

  float *matrix = (float *)memalign(sizeof(float)*4, sizeof(float)*size*size);
  if(matrix==NULL){
    printf("can't allocate the required memory for matrix\n");
    return 0;
  }

  float *matrix_Neo = (float *)memalign(sizeof(float)*4, sizeof(float)*size*size);
  if(matrix_Neo==NULL){
    printf("can't allocate the required memory for matrix\n");
    return 0;
  }

  float *matrix_Trini = (float *)memalign(sizeof(float)*4, sizeof(float)*size*size);
  if(matrix_Trini==NULL){
    printf("can't allocate the required memory for matrix\n");
    return 0;
  }

  float *matrix_Morph = (float *)memalign(sizeof(float)*4, sizeof(float)*size*size);
  if(matrix_Morph==NULL){
    printf("can't allocate the required memory for matrix\n");
    return 0;
  }

  // Init matrices for default init gives not all zero
  init_Matrix(matrix_Neo, size);
  init_Matrix(matrix_Trini, size);
  init_Matrix(matrix_Morph, size);

  matrix_vector_gen(size, matrix_Neo);
  matrix_vector_gen(size, matrix_Trini);

  matrix_mult_sse(size, matrix_Neo, matrix_Trini, matrix_Morph);

  print_Matrix(matrix_Morph, size);

  free(matrix);
  free(matrix_Neo);
  free(matrix_Trini);
  free(matrix_Morph);
  return 1;
}
