#!/usr/bin/bash

iterations=10 # edit also the ../plot_speedup.py script

echo " " > data/speedup.txt

N=$(nproc --all)
N=$(($N + 3))

for i in $(seq 1 $N)
do
    echo "Running with $i threads over $N ($(nproc --all) logical core + 3 hyperthreading)"
    for j in $(seq 1 $iterations)
    do
        echo "   Iteration $j on $iterations"
        export OMP_NUM_THREADS=$i
        { time ./data/main.out ; } 2>> data/speedup.txt
        sleep 2
    done
done