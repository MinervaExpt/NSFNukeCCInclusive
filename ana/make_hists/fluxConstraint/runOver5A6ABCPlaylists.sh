#!/bin/bash

for playlist in minervame5A minervame6A minervame6B minervame6C 
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done