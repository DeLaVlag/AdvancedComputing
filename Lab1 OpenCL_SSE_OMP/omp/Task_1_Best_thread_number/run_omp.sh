#!/bin/sh

thread_num=$1
matrix_size=$2

#echo "compile application"

gcc -g -fopenmp -lm matrix_omp.c -o matrix_omp.exe

#echo "setting up number of threads value"
export OMP_NUM_THREADS=$thread_num

#echo "executing the application"
./matrix_omp.exe $matrix_size

rm -fr *~ matrix_omp.exe
