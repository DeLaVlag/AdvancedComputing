#!/bin/sh

for number in $(seq 1 1 3)
do
	LOGFILE="Lab1PartBTask5-$number.log"
	sh run_bonus.sh > $LOGFILE
done