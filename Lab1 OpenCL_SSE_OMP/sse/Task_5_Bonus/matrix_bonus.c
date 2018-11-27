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
void matrix_vector_gen(int size, float *matrix, float *vector) {
	int i;
	for (i = 0; i < size; i++)
		vector[i] = ((float)rand()) / 65535.0f;
	for (i = 0; i < size * size; i++)
		matrix[i] = ((float)rand()) / 5307.0f;
}

/****************************************************
the following function calculate the below equation
   vector_out = vector_in x matrix_in
 ***************************************************/
void matrix_mult_sq(int size, float *vector_in,
                    float *matrix_in, float *vector_out) {
	int rows, cols;
	int j;
	for (cols = 0; cols < size; cols++) {
		vector_out[cols] = 0.0;
		for (j = 0, rows = 0; rows < size; j++, rows++) {
			vector_out[cols] += vector_in[j] * matrix_in[rows * size + cols];
		}
	}
}

/****************************************************
the following function calculates the below equation
   vector_out = vector_in x matrix_in
   With the SSE parallellization.
 ***************************************************/
void matrix_mult_sse(int size, float *vector_in,
                     float *matrix_in, float *vector_out) {
	__m128 a_line, b_line, r_line;
	int i, j;
	for (i = 0; i < size; i += 4) {
		j = 0;
		b_line = _mm_load_ps(&matrix_in[i]);               // b_line = vec4(matrix[i][0])
		a_line = _mm_set1_ps(vector_in[j]);                // a_line = vec4(vector_in[0])
		r_line = _mm_mul_ps(a_line, b_line);               // r_line = a_line * b_line
		for (j = 1; j < size; j++) {
			b_line = _mm_load_ps(&matrix_in[j * size + i]);// a_line = vec4(column(a, j))
			a_line = _mm_set1_ps(vector_in[j]);            // b_line = vec4(b[i][j])
			r_line = _mm_add_ps(_mm_mul_ps(a_line, b_line), r_line); // r_line += a_line * b_line
		}
		_mm_store_ps(&vector_out[i], r_line);              // r[i] = r_line
	}
}

/****************************************************
the following function calculates the below equation
   vector_out = vector_in x matrix_in
   With the OMP parallellization.
 ***************************************************/
void matrix_mult_omp(int size, float *vector_in,
                     float *matrix_in, float *vector_out) {
	int rows, cols;
	int j;
	# pragma omp parallel       \
	shared(size, vector_in, matrix_in, vector_out)  \
	private(rows, cols, j)
	# pragma omp for
	for (cols = 0; cols < size; cols++) {
		vector_out[cols] = 0.0;
		for (j = 0, rows = 0; rows < size; j++, rows++) {
			vector_out[cols] += vector_in[j] * matrix_in[rows * size + cols];
		}
	}
}

/****************************************************
the following function calculates the below equation
   vector_out = vector_in x matrix_in
   With the OMP and SSE parallellization.
 ***************************************************/
void matrix_mult_both(int size, float *vector_in,
                      float *matrix_in, float *vector_out) {

	__m128 a_line, b_line, r_line;
	int i, j;
	# pragma omp parallel       \
	shared(size, vector_in, matrix_in, vector_out)  \
	private(a_line, b_line, r_line,i,j)
	# pragma omp for
	for (i = 0; i < size; i += 4) {
		j = 0;
		b_line = _mm_load_ps(&matrix_in[i]);               // b_line = vec4(matrix[i][0])
		a_line = _mm_set1_ps(vector_in[j]);                // a_line = vec4(vector_in[0])
		r_line = _mm_mul_ps(a_line, b_line);               // r_line = a_line * b_line
		for (j = 1; j < size; j++) {
			b_line = _mm_load_ps(&matrix_in[j * size + i]);// a_line = vec4(column(a, j))
			a_line = _mm_set1_ps(vector_in[j]);            // b_line = vec4(b[i][j])
			r_line = _mm_add_ps(_mm_mul_ps(a_line, b_line), r_line); // r_line += a_line * b_line
		}
		_mm_store_ps(&vector_out[i], r_line);              // r[i] = r_line
	}
}

/****************************************************
the following function checks if two matrices are
the same.
 ***************************************************/
int matrix_check(int size, float *matrixA, float *matrixB) {
	//check
	int i;
	for (i = 0; i < size; i++) {
		if ((int)matrixA[i] != (int)matrixB[i]) {
			printf("wrong at position %d\n", i);
			return 0;
		}
	}
}

/****************************************************
the main-function
 ***************************************************/
int main(int argc, char *argv[]) {
	if (argc < 2) {
		printf("Usage: %s matrix/vector_size\n", argv[0]);
		return 0;
	}

	int size = atoi(argv[1]);
	if (size % 4 != 0) {
		printf("This version implements for ""size = 4*n"" only\n");
		return 0;
	}

	///Declare aligned memory
	float *vector = (float *)memalign(sizeof(float) * 4, sizeof(float) * size);
	if (vector == NULL) {
		printf("can't allocate the required memory for vector\n");
		return 0;
	}

	float *matrix = (float *)memalign(sizeof(float) * 4, sizeof(float) * size * size);
	if (matrix == NULL) {
		printf("can't allocate the required memory for matrix\n");
		free(vector);
		return 0;
	}

	float *result_sq = (float *)memalign(sizeof(float) * 4, sizeof(float) * size);
	if (result_sq == NULL) {
		printf("can't allocate the required memory for result_sq\n");
		free(vector);
		free(matrix);
		return 0;
	}

	float *result_omp = (float *)memalign(sizeof(float) * 4, sizeof(float) * size);
	if (result_omp == NULL) {
		printf("can't allocate the required memory for result_pl\n");
		free(vector);
		free(matrix);
		free(result_omp);
		free(result_sq);
		return 0;
	}

	float *result_sse = (float *)memalign(sizeof(float) * 4, sizeof(float) * size);
	if (result_sse == NULL) {
		printf("can't allocate the required memory for result_pl\n");
		free(vector);
		free(matrix);
		free(result_omp);
		free(result_sq);
		free(result_sse);
		return 0;
	}

	float *result_both = (float *)memalign(sizeof(float) * 4, sizeof(float) * size);
	if (result_both == NULL) {
		printf("can't allocate the required memory for result_pl\n");
		free(vector);
		free(matrix);
		free(result_omp);
		free(result_sq);
		free(result_sse);
		free(result_both);
		return 0;
	}

	///Generate a vector and matrix with random values
	matrix_vector_gen(size, matrix, vector);

	double time_sq;
	double time_sse;
	double time_omp;
	double time_both;

	///Get time performance of normal exectuion
	time_sq = omp_get_wtime();
	matrix_mult_sq(size, vector, matrix, result_sq);
	time_sq = omp_get_wtime() - time_sq;

	///Get time performance of SSE exectuion
	time_omp = omp_get_wtime();
	matrix_mult_omp(size, vector, matrix, result_omp);
	time_omp = omp_get_wtime() - time_omp;
	matrix_check(size, result_omp, result_sq);

	///Get time performance of SSE exectuion
	time_sse = omp_get_wtime();
	matrix_mult_sse(size, vector, matrix, result_sse);
	time_sse = omp_get_wtime() - time_sse;
	matrix_check(size, result_sse, result_sq);

	///Get time performance of both OMP combined with SSE
	time_both = omp_get_wtime();
	matrix_mult_both(size, vector, matrix, result_both);
	time_both = omp_get_wtime() - time_both;
	matrix_check(size, result_both, result_sq);

	printf("SZ:\t%04d\t",size);
	printf("SEQ:\t%f\t", time_sq);
	printf("OMP:\t%f\t", time_omp);
	printf("SSE:\t%f\t", time_sse);
	printf("BOTH:\t%f\t", time_both);
	printf("SEQ / PAR:\t%f\n", time_sq / time_both);

	free(vector);
	free(matrix);
	free(result_omp);
	free(result_sq);
	free(result_sse);
	free(result_both);
	return 1;
}
