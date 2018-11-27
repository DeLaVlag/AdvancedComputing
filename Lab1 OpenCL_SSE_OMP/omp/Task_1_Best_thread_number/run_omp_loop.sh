#!/bin/sh

for matrix_size in $(seq 2000 2000 3000)
do
	for thread_num in $(seq 1 1 20)
	do
		echo ""
		echo "$matrix_size"
		echo "$thread_num"
	#	echo -n "S:\t\t"
	#	./run_sse.sh $matrix_size
	#	echo -n "P:\t\t"
		for numb in $(seq 1 1 15)
		do
			./run_omp.sh $thread_num $matrix_size
			echo ""
		done
	done
	echo ""
done

