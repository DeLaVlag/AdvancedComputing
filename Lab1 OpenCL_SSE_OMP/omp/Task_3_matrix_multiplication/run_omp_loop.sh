#!/bin/sh

numRepeat=5

for matrix_size in $(seq 1000 1000 5000)
do
	#compare matrix vector multiplication to matrix matrix multiplation

	echo  "---------------------------------------------\n"
	echo  " Matrix vector Multiplication "
	echo  "---------------------------------------------"

	for thread_num in $(seq 8 2 9)
	do
		echo -n "MS: $matrix_size "
		echo -n "TN: $thread_num "
		echo  ""
		for numb in $(seq 1 1 $numRepeat)
		do
			./run_omp_vTimesm.sh $thread_num $matrix_size
			echo ""
		done
	done
	echo ""

	echo  "---------------------------------------------"
	echo  " Matrix Multiplication "
	echo  "---------------------------------------------"

	for thread_num in $(seq 8 2 9)
	do
		echo -n "MS: $matrix_size "
		echo -n "TN: $thread_num "
		echo  ""
		for numb in $(seq 1 1 $numRepeat)
		do
			./run_omp_mTimesm.sh $thread_num $matrix_size
			echo ""
		done
	done
	echo ""

done

