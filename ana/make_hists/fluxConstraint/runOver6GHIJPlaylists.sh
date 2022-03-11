#!/bin/bash

for playlist in minervame6G minervame6H minervame6I minervame6J 
#for playlist in minervame6G minervame6H minervame6J
do
    echo $playlist
    ./runEventLoop_trueCCwo2p2h_optimised_petal . $playlist
done