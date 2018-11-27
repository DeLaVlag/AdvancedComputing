#!/bin/sh

rows_size=$1
cols_size=$2


# echo "compile application"

gcc -g -lm -msse -fopenmp  matrix_sse.c -o matrix_sse.exe

# echo "executing the application"
./matrix_sse.exe $rows_size $cols_size

rm -fr *~ matrix_sse.exe
