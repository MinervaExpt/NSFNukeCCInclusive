#!/bin/bash

for playlist in minervame6D minervame6E minervame6F 
#for playlist in minervame6F 
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done