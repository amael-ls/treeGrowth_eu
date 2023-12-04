#!/bin/bash

echo "You are running $0"
echo "Make sure that there are 6 CPUs available."
echo "You have 10 seconds to cancel $0 bash script in case you did not check CPUs"

sleep 10

ls_species=("Betula pendula" "Fagus sylvatica" "Picea abies" "Pinus pinaster" "Pinus sylvestris" "Quercus petraea")

for species_id in "${ls_species[@]}"
do
	echo ${species_id}
	/usr/bin/Rscript --vanilla get_optima.R "${species_id}" 1 > optima_"${species_id}".txt 2>&1 &
done
