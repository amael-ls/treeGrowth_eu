#!/bin/bash

echo "You are running $0. If you have not done it yet, maybe you should check how many species you are trying to run"
echo "Make sure that there are enough CPUSs left available for others"
echo "You have 10 seconds to cancel $0 bash script in case you did not check CPUs or would like to change the species runned"

sleep 10

echo "Starting first run, the other runs will start 45 seconds later (to let enough time to the first run to create the folders)"

today=$(date +"%Y-%m-%d_%T")

for species_id in {1..45} # Maybe run by chunk of 5, except if your computer can handle 45 species x 4 runs x 4 chains!
do
	/usr/bin/Rscript --vanilla growth_subsample.R $species_id 1 15000 > ${species_id}_1_${today}.txt 2>&1 &
done

sleep 45

for species_id in {1..45} # Maybe run by chunk of 5, except if your computer can handle 45 species x 4 runs x 4 chains!
do
	for run_id in {2..4}
	do
		/usr/bin/Rscript --vanilla growth_subsample.R $species_id $run_id 15000 > ${species_id}_${run_id}_${today}.txt 2>&1 &
	done
done
