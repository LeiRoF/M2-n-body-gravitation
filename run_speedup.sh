#!/usr/bin/bash

echo "⚙️ Compiling..."
./src/compile.sh
echo "✅ Done!"

echo " " > data/tmp_speedup.txt

N=$(nproc --all)
N=$(($N + 3))

echo "🏃 Running simulation with 1 to  $N threads ($(nproc --all) logical core + 3 hyperthreading)"
for i in $(seq 1 $N)
do
    echo "   Running with $i threads"
    export OMP_NUM_THREADS=$i
    { time ./data/main.out ; } 2>> data/tmp_speedup.txt
    sleep 2
done
echo "✅ Done!"

python3 ./src/format_speedup.py

rm data/tmp_speedup.txt