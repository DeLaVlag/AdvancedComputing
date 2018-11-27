#!/bin/sh

for matrix_size in $(seq 1000 100 10000)
do
#	echo ""
#	echo matrix_size
	printf "$matrix_size\t"
	./run_sse.sh $matrix_size
#	echo ""
#	echo ""
done

