#!/bin/bash

for playlist in minervame1A minervame1B minervame1C minervame1D minervame1E minervame1F minervame1G minervame1L minervame1M minervame1N minervame1O minervame1P
do
    echo $playlist
    ./runEventLoopALL_trueCCwo2p2h . 3 26 $playlist
done