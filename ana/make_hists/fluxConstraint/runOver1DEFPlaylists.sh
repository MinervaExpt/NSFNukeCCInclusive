#!/bin/bash

#for playlist in minervame1D minervame1E minervame1F
for playlist in minervame1F
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done