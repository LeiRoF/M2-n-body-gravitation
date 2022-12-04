export OMP_NUM_THREADS=$(nproc --all)
gfortran -O3 -fopenmp -Wall -fbounds-check ./main.f90 -o ./data/main.out