#!/bin/bash

for playlist in minervame1G minervame1L minervame1M
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised . $playlist
done