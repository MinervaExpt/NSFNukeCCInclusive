#!/bin/bash

#for playlist in minervame1A minervame1B minervame1C minervame1D minervame1E minervame1F minervame1G minervame1L minervame1M minervame1N minervame1O minervame1P
#for playlist in minervame5A minervame6A minervame6B minervame6C minervame6D minervame6E minervame6F minervame6I minervame6G minervame6H minervame6I minervame6J
for playlist in minervame1E
do 
    echo $playlist
    find "${PWD}" /pnfs/minerva/persistent/MAD_ntuples/production_p1/MuonKuldged/Neutrinos/FullDetector/Merged_mc_ana_me1E_FirstProductionAnatuples_DualVertex/ -name 'MasterAnaDev_mc_AnaTuple_run*.root' > ${playlist}_mc_DualVertex_FullDetector.txt
    #xrdify $playlist.txt > MasterAnaDev_MC_${playlist}_MuonKludged.txt
done