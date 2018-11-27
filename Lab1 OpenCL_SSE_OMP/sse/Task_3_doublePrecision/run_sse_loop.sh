#!/bin/sh

for rows_size in $(seq 1000 100 1000000)
do
#	echo ""
#	echo matrix_size
	printf "$rows_size\t"
	./run_sse.sh $rows_size $rows_size
#	echo ""
#	echo ""
done

