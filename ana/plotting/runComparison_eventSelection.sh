#!/bin/bash
echo tune comparison + chi2/ndf for bkg subtracted distributions

echo Carbon
python plotCVError_dataMCratio_comparison.py 3 06 1

echo Iron
python plotCVError_dataMCratio_comparison.py 2 26 1
python plotCVError_dataMCratio_comparison.py3 26 1
python plotCVError_dataMCratio_comparison.py 5 26 1

echo Lead
python plotCVError_dataMCratio_comparison.py 2 82 1
python plotCVError_dataMCratio_comparison.py 3 82 1
python plotCVError_dataMCratio_comparison.py 4 82 1
python plotCVError_dataMCratio_comparison.py 5 82 1

echo Tracker
python plotCVError_dataMCratio_comparison.py 99 99 1