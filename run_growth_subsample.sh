#!/bin/bash

echo "You are running $0. If you have not done it yet, maybe use tmux to open a new window"
echo "Make sure that there are at least 30 CPUs available."
echo "You have 10 seconds to cancel $0 bash script in case you did not check CPUs or would like to use tmux"

sleep 10

echo "Starting first run, the other runs will start 45 seconds later in background"

today=$(date +"%Y-%m-%d_%T")

for species_id in 2 5 17 48
do
	/usr/bin/Rscript --vanilla growth_subsample.R $species_id 1 15000 > ${species_id}_1_${today}.txt 2>&1 &
done

sleep 45

for species_id in 17
do
	for run_id in {2..4}
	do
		/usr/bin/Rscript --vanilla growth_subsample.R $species_id $run_id 15000 > ${species_id}_${run_id}_${today}.txt 2>&1 &
	done
done
