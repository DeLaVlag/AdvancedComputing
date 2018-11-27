#!/bin/sh


	for matrix_size in $(seq 1000 100 10000) #loop matrix size
	do
		for thread_num in $(seq 8 1 8) #thread loop
		do
			echo ""
			echo "$matrix_size"
			echo "$thread_num"
			
			for numb in $(seq 1 1 1)
			do
				./run_omp.sh $thread_num $matrix_size
				echo ""
			done
		done
		echo ""
	done


