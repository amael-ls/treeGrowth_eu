#!/bin/bash

echo "You are running $0. If you have not done it yet, maybe use tmux to open a new window"
echo "Make sure that there are 48 CPUs available."
echo "You have 10 seconds to cancel $0 bash script in case you did not check CPUs or would like to use tmux"

sleep 10

today=$(date +"%Y-%m-%d_%T")

echo "Running at $today"

ls_species=("Betula pendula" "Fagus sylvatica" "Picea abies" "Pinus pinaster" "Pinus sylvestris" "Quercus petraea")

for species_id in "${ls_species[@]}"
do
	echo ${species_id}
	/usr/bin/Rscript --vanilla waic_sa.R $species_id 1 min_max > "${species_id}"_${today}.txt 2>&1 &
done
