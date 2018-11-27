#!/bin/sh

for number in $(seq 1 1 4)
do
	LOGFILE="Lab1PartATask2-$number.log"
	sh run_omp_loop.sh > $LOGFILE
done