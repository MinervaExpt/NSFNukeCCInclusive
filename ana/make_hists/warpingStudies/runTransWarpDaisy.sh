#!/bin/bash

#Usage: runTransWarp.sh 
PWD=$1 #infile dir
origPWD=$2 #infile dir
out=$3 # outfile dir
targetID=$4 # targetID
targetZ=$5 # target Z
plistTag=$6 # What tag
warp=$7

# h_bkg_subtracted_mc_daisy_9_x

for VARIABLE in Enu x #pTmu pZmu
do
    echo For petal ${VARIABLE} start:
    for PETAL in 0 1 2 3 4 5 6 7 8 9 10 11
    do
        echo For petal ${PETAL} start:
        OUTFILE_NAME=${warp}_${VARIABLE}_${PETAL}

        FAKE_DATA_RECO_FILE=${PWD}/BackgroundSubtracted/BkgSubtracted_EventSelection_daisy_${plistTag}_t${targetID}_z${targetZ}_sys.root
        FAKE_DATA_RECO_HISTO=h_bkg_subtracted_mc_daisy_${PETAL}_${VARIABLE}
        # Must be selected signal!! - background substracted
        # for tracker, this is the selected signal in event selection
        # for targets, this is slightly different due to plastic sideband tuning

        FAKE_DATA_TRUTH_FILE=${PWD}/Efficiency/Efficiency_Daisy_${plistTag}_t${targetID}_z${targetZ}_sys.root  
        FAKE_DATA_TRUTH_HISTO=h_mc_daisy_${PETAL}_${VARIABLE}

        ORIG_RECO_FILE=${origPWD}/BackgroundSubtracted/BkgSubtracted_EventSelection_daisy_${plistTag}_t${targetID}_z${targetZ}_sys.root
        ORIG_RECO_HISTO=h_bkg_subtracted_mc_daisy_${PETAL}_${VARIABLE}

        ORIG_TRUTH_FILE=${origPWD}/Efficiency/Efficiency_Daisy_${plistTag}_t${targetID}_z${targetZ}_sys.root
        ORIG_TRUTH_HISTO=h_mc_daisy_${PETAL}_${VARIABLE}

        MIGRATION_FILE=${origPWD}/Migration/Migration_Daisy_${plistTag}_t${targetID}_z${targetZ}_sys.root
        MIGRATION_HISTO=selected_Migration_daisy_${PETAL}_${VARIABLE}
        #selected_Migration_${VARIABLE}


        TransWarpExtraction --output_file ${out}/Warping_t${targetID}_z${targetZ}_${plistTag}_$OUTFILE_NAME.root --data_file $FAKE_DATA_RECO_FILE --data $FAKE_DATA_RECO_HISTO --data_truth_file $FAKE_DATA_TRUTH_FILE  --data_truth $FAKE_DATA_TRUTH_HISTO  --reco_file $ORIG_RECO_FILE --reco $ORIG_RECO_HISTO --truth_file $ORIG_TRUTH_FILE --truth $ORIG_TRUTH_HISTO --migration_file $MIGRATION_FILE --migration $MIGRATION_HISTO --num_iter 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 -C 2 --max_chi2 200 --num_uni 100

    done
    echo For variable ${VARIABLE} DONE  
done
