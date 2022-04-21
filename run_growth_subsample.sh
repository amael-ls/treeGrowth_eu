#!/bin/bash

echo "You are running $0. If you have not done it yet, maybe use tmux to open a new window"
echo "Make sure that there are at least 30 CPUs available."
echo "You have 20 seconds to cancel $0 bash script in case you did not check CPUs or would like to use tmux"

sleep 20

echo "Starting first run, the other runs will start 45 seconds later in background"

for species_id in 17 48
do
	Rscript --vanilla growth_subsample.R $species_id 1 8000 &
done

sleep 45

for species_id in 17 48
do
	for run_id in {2..4}
	do
		Rscript --vanilla growth_subsample.R $species_id $run_id 8000 &
	done
done
