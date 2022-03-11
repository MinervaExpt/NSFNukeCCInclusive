#!/bin/bash

for playlist in minervame1A minervame1B minervame1C
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done