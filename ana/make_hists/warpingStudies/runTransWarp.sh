#!/bin/bash

#Usage: runTransWarp.sh 

VARIABLE=x
OUTFILE_NAME=Orig_x

PWD=/minerva/data2/users/anezkak/TrackerWarpingStudies/DeuteriumGeniePiTune/

FAKE_DATA_RECO_FILE=${PWD}Hists_EventSelectionTracker_ML_ME6A_nosys_t99_z99_AntiNu.root
FAKE_DATA_RECO_HISTO=selected_mc_reco_signal_${VARIABLE}
# Must be selected signal!!

FAKE_DATA_TRUTH_FILE=${PWD}Hists_Efficiency_ML_ME6A_nosys_t99_z99_AntiNu.root
FAKE_DATA_TRUTH_HISTO=h_mc_${VARIABLE}

ORIG_RECO_FILE=${PWD}Hists_EventSelectionTracker_ML_ME6A_nosys_t99_z99_AntiNu.root
ORIG_RECO_HISTO=selected_mc_reco_signal_${VARIABLE}

ORIG_TRUTH_FILE=${PWD}Hists_Efficiency_ML_ME6A_nosys_t99_z99_AntiNu.root
ORIG_TRUTH_HISTO=h_mc_${VARIABLE}

MIGRATION_FILE=${PWD}Hists_Migration_ML_ME6A_nosys_t99_z99_AntiNu.root
MIGRATION_HISTO=response1d_${VARIABLE}_migration
#selected_Migration_${VARIABLE}


TransWarpExtraction --output_file Warping_Tracker_NukeCC_ME6A_$OUTFILE_NAME.root --data_file $FAKE_DATA_RECO_FILE --data $FAKE_DATA_RECO_HISTO --data_truth_file $FAKE_DATA_TRUTH_FILE  --data_truth $FAKE_DATA_TRUTH_HISTO  --reco_file $ORIG_RECO_FILE --reco $ORIG_RECO_HISTO --truth_file $ORIG_TRUTH_FILE --truth $ORIG_TRUTH_HISTO --migration_file $MIGRATION_FILE --migration $MIGRATION_HISTO --num_iter 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 -C 2 --max_chi2 200 --num_uni 100
