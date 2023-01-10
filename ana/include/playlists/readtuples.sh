#!/bin/bash

#for playlist in minervame1A minervame1B minervame1C minervame1D minervame1E minervame1F minervame1G minervame1L minervame1M minervame1N minervame1O minervame1P
#for playlist in minervame5A minervame6A minervame6B minervame6C minervame6D minervame6E minervame6F minervame6I minervame6G minervame6H minervame6I minervame6J
for playlist in minervame6A
do 
    echo $playlist
    # MC
    find "${PWD}" /pnfs/minerva/persistent/users/zdar/Merged_mc_ana_me6A_DualVertex_p3/ -name 'MasterAnaDev_mc_AnaTuple_run*.root' > ${playlist}_mc_DualVertex_FullDetector_p3.txt
    xrdify ${playlist}_mc_DualVertex_FullDetector_p3.txt > MasterAnaDev_MC_${playlist}.txt
    # Data
    find "${PWD}" /pnfs/minerva/persistent/users/zdar/Merged_data_ana_me6A_DualVertex_p3/ -name 'MasterAnaDev_data_AnaTuple_run*.root' > ${playlist}_data_DualVertex_p3.txt
    xrdify ${playlist}_data_DualVertex_p3.txt > MasterAnaDev_Data_${playlist}.txt

done