#!/bin/bash

for playlist in minervame1N minervame1O minervame1P
do
    echo $playlist
    ./runEventLoopALL_trueCCwo2p2h_tracker . $playlist
done