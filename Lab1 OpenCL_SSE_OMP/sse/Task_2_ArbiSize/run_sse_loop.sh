#!/bin/sh

for rows_size in $(seq 1000 101 10000)
do
	# cols_size = rows_size
#	echo ""
#	echo matrix_size
	printf "$rows_size\t"
	./run_sse.sh $rows_size $rows_size
#	echo ""
#	echo ""
done

