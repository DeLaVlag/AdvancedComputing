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


#define PRINTMATRIX(X,S) //print_Matrix((X), (S))

/*****************************************************
the following function generates a "size"-element vector
and a "size x size" matrix
 ****************************************************/
void matrix_vector_gen(int size, double *matrixA, double *matrixB) {
	int i;
	///Generate two matrices with random values
	for (i = 0; i < size * size; i++)
		matrixA[i] = ((double)rand()) / 5307.0;
	for (i = 0; i < size * size; i++)
		matrixB[i] = ((double)rand()) / 5307.0;
}

/****************************************************
the following function calculate the below equation
   matrix_out = matrix_A x matrix_B
 ***************************************************/
void matrix_mult_sq(int size, double *matrixA,
                    double *matrixB, double *matrix_out) {
	int rows1, rows2, cols;
	for (rows1 = 0; rows1 < size; rows1++) {
		for (cols = 0; cols < size; cols++) {
			matrix_out[rows1 * size + cols] = 0.0;
			for (rows2 = 0; rows2 < size; rows2++) {
				matrix_out[rows1 * size + cols] += matrixA[rows1 * size + rows2] * matrixB[rows2 * size + cols];
			}
		}
	}
}

/****************************************************
the following function calculates the equation below
   matrix_out = matrix_A x matrix_B
   With the OMP parallellization.
 ***************************************************/
void matrix_mult_pl(int size, double *matrixA,
                    double *matrixB, double *matrix_out) {
	int rows1, rows2, cols;
	/* shared() the input parameters and output variabel*/
	/* privat() the local variables */
	# pragma omp parallel       \
	shared(size, matrixA, matrixB, matrix_out)  \
	private(rows1, rows2, cols)
	# pragma omp for
	for (rows1 = 0; rows1 < size; rows1++) {
		for (cols = 0; cols < size; cols++) {
			matrix_out[rows1 * size + cols] = 0.0;
			for (rows2 = 0; rows2 < size; rows2++) {
				matrix_out[rows1 * size + cols] += matrixA[rows1 * size + rows2] * matrixB[rows2 * size + cols];
			}
		}
	}
}


/****************************************************
function prints a Matrix

 ***************************************************/
void print_Matrix(double *matrix_in, int size) {
	int i = 0, j;

	printf("\n");
	for (j = 0; j < size * size; j++) {
		printf(" %1.3f ", matrix_in[j]);
		i++;
		if (i == size) {
			printf("\n");
			i = 0;
		}
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

	///Create memory for the matrix calculations
	int size = atoi(argv[1]);
	double *matrixA   = (double *)malloc(sizeof(double) * size * size);
	double *matrixB   = (double *)malloc(sizeof(double) * size * size);
	double *result_sq = (double *)malloc(sizeof(double) * size * size);
	double *result_pl = (double *)malloc(sizeof(double) * size * size);
	matrix_vector_gen(size, matrixA, matrixB);
	PRINTMATRIX(matrixB, size);

	double time_sq = 0;
	double time_pl = 0;

	time_sq = omp_get_wtime();
	matrix_mult_sq(size, matrixA, matrixB, result_sq);
	time_sq = omp_get_wtime() - time_sq;
	PRINTMATRIX(result_sq, size);

	time_pl = omp_get_wtime();
	matrix_mult_pl(size, matrixA, matrixB, result_pl);
	time_pl = omp_get_wtime() - time_pl;
	PRINTMATRIX(result_pl, size);

	///Calculate the performance differences between normal execution and
	///This section does not influence the performance measurement.
	printf("SEQ: \t  %f (sec)\t", time_sq);
	printf("PAR: \t(t:%d) \t (p:%d) \t %f (sec)\t",
	       omp_get_max_threads(), omp_get_num_procs(), time_pl);
	printf("SEQ / PAR: \t  %f ", time_sq / time_pl);

	//check, the calculation SER. vs. PAR.
	int i;
	for (i = 0; i < size; i++) {
		if (result_sq[i] != result_pl[i]) {
			printf("wrong at position %d\n", i);
			return 0;
		}
	}

	free(matrixA);
	free(matrixB);
	free(result_sq);
	free(result_pl);
	return 1;
}
