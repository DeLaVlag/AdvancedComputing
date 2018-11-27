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
void matrix_vector_gen(int size, float *matrix, float *vector){
  int i;
  for(i=0; i<size; i++)
    vector[i] = i*1.2f + 1;//((float)rand())/65535.0f;
  for(i=0; i<size*size; i++)
    matrix[i] = i*1.3f + 1;//((float)rand())/5307.0f;
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

void matrix_mult_sse(int size, float *vector_in,
		      float *matrix_in, float *vector_out){
  __m128 a_line, b_line, r_line;
  int i, j;
  for (i=0; i<size; i+=4){
    j = 0;
    b_line = _mm_load_ps(&matrix_in[i]); // b_line = vec4(matrix[i][0])
    a_line = _mm_set1_ps(vector_in[j]);      // a_line = vec4(vector_in[0])
    r_line = _mm_mul_ps(a_line, b_line); // r_line = a_line * b_line
    for (j=1; j<size; j++) {
      b_line = _mm_load_ps(&matrix_in[j*size+i]); // a_line = vec4(column(a, j))
      a_line = _mm_set1_ps(vector_in[j]);  // b_line = vec4(b[i][j])
                                     // r_line += a_line * b_line
      r_line = _mm_add_ps(_mm_mul_ps(a_line, b_line), r_line);
    }
    _mm_store_ps(&vector_out[i], r_line);     // r[i] = r_line
  }
}


int main(int argc, char *argv[]){
  if(argc < 2){
    printf("Usage: %s matrix/vector_size\n", argv[0]);
    return 0;
  }

  int size = atoi(argv[1]);
  if(size%4 != 0){
    printf("This version implements for ""size = 4*n"" only\n");
    return 0;
  }

  float *vector = (float *)memalign(sizeof(float)*4, sizeof(float)*size);//(float *)malloc(sizeof(float)*size);
  if(vector==NULL){
    printf("can't allocate the required memory for vector\n");
    return 0;
  }

  float *matrix = (float *)memalign(sizeof(float)*4, sizeof(float)*size*size);
  if(matrix==NULL){
    printf("can't allocate the required memory for matrix\n");
    free(vector);
    return 0;
  }

  float *result_sq = (float *)memalign(sizeof(float)*4, sizeof(float)*size);
  if(result_sq==NULL){
    printf("can't allocate the required memory for result_sq\n");
    free(vector);
    free(matrix);
    return 0;
  }

  float *result_pl = (float *)memalign(sizeof(float)*4, sizeof(float)*size);
  if(result_pl==NULL){
    printf("can't allocate the required memory for result_pl\n");
    free(vector);
    free(matrix);
    free(result_sq);
    return 0;
  }

  matrix_vector_gen(size, matrix, vector);

  double time_sq;
  double time_sse;

  time_sq = omp_get_wtime();
  matrix_mult_sq(size, vector, matrix, result_sq);
  time_sq = omp_get_wtime() - time_sq;

  time_sse = omp_get_wtime();
  matrix_mult_sse(size, vector, matrix, result_pl);
  time_sse = omp_get_wtime() - time_sse;

  printf("%f:", time_sq);
  printf("%f:\n", time_sse);
  //printf("SEQ: %f (sec)\t", time_sq);
  //printf("PAR: %f (sec)\t", time_sse);
  //printf("SEQ / PAR: %f \n", time_sq/time_sse);

  //printf("SEQUENTIAL EXECUTION: %f (sec)\n",time_sq);
  //printf("PARALLEL EXECUTION: %f (sec)\n", time_sse);

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
  free(matrix);
  free(result_sq);
  free(result_pl);
  return 1;
}
