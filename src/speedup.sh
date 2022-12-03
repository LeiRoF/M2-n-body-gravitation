#!/bin/bash

echo " " > data/speedup.txt

for i in $(seq 1 $(nproc --all))
do
    for j in {1..10}
    do
        export OMP_NUM_THREADS=$i
        { time ./data/main.out ; } 2>> data/speedup.txt
    done
done