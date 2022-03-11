#!/bin/bash

#minerva/data2/users/anezkak/flux4Daisy_files_2022/validation2

STR="/minerva/data2/users/anezkak/fluxFinal/Antineutrinos"
echo $STR


python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_trackerC_Enu tracker |  tee $STR/tracker_flux.txt
python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_t25pb_Enu lead | tee $STR/lead25_flux.txt
python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_t15pb_Enu lead_t15 | tee $STR/lead15_flux.txt
python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_t3c_Enu carbon | tee $STR/carbon_flux.txt
python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_t25fe_Enu iron | tee $STR/iron_flux.txt
python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_t15fe_Enu iron_t15 | tee $STR/iron15_flux.txt
python flux_antinu_calculation.py $STR $STR/CombinedPlaylists_AntiNu_sys_Enu_vs_Petal.root selected_mc_truth_waterO_Enu water_apothem | tee $STR/water_apothem_flux.txt
