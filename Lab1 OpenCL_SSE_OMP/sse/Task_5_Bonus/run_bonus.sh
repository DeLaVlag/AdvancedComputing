#!/bin/sh

thread_num=8

for matrix_size in $(seq 0 100 10000)
do

	gcc -g -lm -msse -fopenmp  matrix_bonus.c -o matrix_bonus.exe

	export OMP_NUM_THREADS=$thread_num

	./matrix_bonus.exe $matrix_size 

	rm -fr *~ matrix_bonus.exe

done