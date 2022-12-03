export OMP_NUM_THREADS=$(nproc --all)

gfortran -O3 -fopenmp ./config.f90 ./src/main.f90 -o ./data/main.out

rm ./config.mod
# rm ./src/settings.mod