#!/bin/bash

for playlist in minervame1N minervame1O minervame1P
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done