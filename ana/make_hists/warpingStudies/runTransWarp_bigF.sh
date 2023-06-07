#!/bin/bash

#Usage: runTransWarp.sh 
PWD=$1 #infile dir
origPWD=$2 #infile dir
out=$3 # outfile dir
targetID=$4 # targetID
targetZ=$5 # target Z
plistTag=$6 # What tag
warp=$7
bigF=$8


for VARIABLE in Enu x pTmu1D pZmu1D ThetamuDeg
do

    #VARIABLE=Enu
    OUTFILE_NAME=${warp}_${VARIABLE}

    FAKE_DATA_RECO_FILE=${PWD}/BackgroundSubtracted/BkgSubtracted_EventSelection_${plistTag}_t${targetID}_z${targetZ}_sys.root
    FAKE_DATA_RECO_HISTO=h_bkg_subtracted_mc_${VARIABLE}
    # Must be selected signal!! - background substracted
    # for tracker, this is the selected signal in event selection
    # for targets, this is slightly different due to plastic sideband tuning

    FAKE_DATA_TRUTH_FILE=${PWD}/Efficiency/Efficiency_${plistTag}_t${targetID}_z${targetZ}_sys.root  
    FAKE_DATA_TRUTH_HISTO=h_mc_${VARIABLE}

    ORIG_RECO_FILE=${origPWD}/BackgroundSubtracted/BkgSubtracted_EventSelection_${plistTag}_t${targetID}_z${targetZ}_sys.root
    ORIG_RECO_HISTO=h_bkg_subtracted_mc_${VARIABLE}

    ORIG_TRUTH_FILE=${origPWD}/Efficiency/Efficiency_${plistTag}_t${targetID}_z${targetZ}_sys.root
    ORIG_TRUTH_HISTO=h_mc_${VARIABLE}

    MIGRATION_FILE=${origPWD}/Migration/Migration_${plistTag}_t${targetID}_z${targetZ}_sys.root
    MIGRATION_HISTO=selected_Migration_${VARIABLE}
    #response1d_${VARIABLE}_migration
    #selected_Migration_${VARIABLE}


    TransWarpExtraction --output_file ${out}/Warping_t${targetID}_z${targetZ}_${plistTag}_$OUTFILE_NAME.root --data_file $FAKE_DATA_RECO_FILE --data $FAKE_DATA_RECO_HISTO --data_truth_file $FAKE_DATA_TRUTH_FILE  --data_truth $FAKE_DATA_TRUTH_HISTO  --reco_file $ORIG_RECO_FILE --reco $ORIG_RECO_HISTO --truth_file $ORIG_TRUTH_FILE --truth $ORIG_TRUTH_HISTO --migration_file $MIGRATION_FILE --migration $MIGRATION_HISTO --num_iter 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 -C 2 --max_chi2 200 --num_uni 100 -F ${bigF}
done
