#!/bin/sh

matrix_size=$1

# echo "compile application"

gcc -g -lm -msse -fopenmp  matrix_sse.c -o matrix_sse.exe

# echo "executing the application"
./matrix_sse.exe $matrix_size

rm -fr *~ matrix_sse.exe
