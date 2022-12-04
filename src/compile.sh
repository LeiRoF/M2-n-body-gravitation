export OMP_NUM_THREADS=$(nproc --all)
gfortran -O3 -fopenmp -Wall -fcheck=bounds ./main.f90 -o ./data/main.out