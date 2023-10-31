#!/bin/bash

echo "You are running $0. If you have not done it yet, maybe run htop"
echo "Make sure that there are at least 30 CPUs available and 400 GB RAM."
echo "You have 10 seconds to cancel $0 bash script in case you did not check CPUs or would like to use tmux"

# sleep 10

echo "Starting"

today=$(date +"%Y-%m-%d_%T")

IFS=""

ls_species=("Betula pendula" "Fagus sylvatica" "Picea abies" "Pinus pinaster" "Pinus sylvestris" "Quercus petraea")

for species in ${ls_species[@]}
do
	echo "Running for $species" 
	/usr/bin/Rscript --vanilla plotGrowth.R $species 1 pr tas > ${species}_1.txt 2>&1 &
done
