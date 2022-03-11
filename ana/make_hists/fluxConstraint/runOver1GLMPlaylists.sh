#!/bin/bash

#for playlist in minervame1G minervame1L minervame1M
for playlist in minervame1M

do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done