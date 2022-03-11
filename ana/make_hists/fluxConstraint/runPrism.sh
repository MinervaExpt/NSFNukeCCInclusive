#!/bin/bash

#minerva/data2/users/anezkak/flux4Daisy_files_2022/validation2
STR="/minerva/data2/users/anezkak/fluxFinal/Neutrinos"
echo $STR

for material in carbon #iron iron_t15 lead lead_t15 water_apothem
do
    echo material
    python flux_prism_2d.py $STR $material 140 | tee $STR/logs/${material}_fit_log.txt
done
