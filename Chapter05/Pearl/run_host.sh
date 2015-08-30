#!/bin/sh
if env | grep -q ^MyNumberOfCores=
then
unset OMP_NUM_THREADS
unset KMP_AFFINITY
unset KMP_PLACE_THREADS
./diffusion_base_host
./diffusion_base_host
./diffusion_base_host

./diffusion_base_Large_host
./diffusion_base_Large_host
./diffusion_base_Large_host

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=$MyNumberOfCores
./diffusion_omp_host
./diffusion_omp_host
./diffusion_omp_host

./diffusion_omp_Large_host
./diffusion_omp_Large_host
./diffusion_omp_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_omp_host
./diffusion_omp_host
./diffusion_omp_host

./diffusion_omp_Large_host
./diffusion_omp_Large_host
./diffusion_omp_Large_host

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=$MyNumberOfCores
./diffusion_ompvect_host
./diffusion_ompvect_host
./diffusion_ompvect_host

./diffusion_ompvect_Large_host
./diffusion_ompvect_Large_host
./diffusion_ompvect_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_ompvect_host
./diffusion_ompvect_host
./diffusion_ompvect_host

./diffusion_ompvect_Large_host
./diffusion_ompvect_Large_host
./diffusion_ompvect_Large_host

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=$MyNumberOfCores
./diffusion_peel_host
./diffusion_peel_host
./diffusion_peel_host

./diffusion_peel_Large_host
./diffusion_peel_Large_host
./diffusion_peel_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_peel_host
./diffusion_peel_host
./diffusion_peel_host

./diffusion_peel_Large_host
./diffusion_peel_Large_host
./diffusion_peel_Large_host

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=$MyNumberOfCores
./diffusion_tiled_host
./diffusion_tiled_host
./diffusion_tiled_host

./diffusion_tiled_Large_host
./diffusion_tiled_Large_host
./diffusion_tiled_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_tiled_host
./diffusion_tiled_host
./diffusion_tiled_host

./diffusion_tiled_Large_host
./diffusion_tiled_Large_host
./diffusion_tiled_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_tiled_HT1_host
./diffusion_tiled_HT1_host
./diffusion_tiled_HT1_host

./diffusion_tiled_HT1_Large_host
./diffusion_tiled_HT1_Large_host
./diffusion_tiled_HT1_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_tiled_HT2_host
./diffusion_tiled_HT2_host
./diffusion_tiled_HT2_host

./diffusion_tiled_HT2_Large_host
./diffusion_tiled_HT2_Large_host
./diffusion_tiled_HT2_Large_host

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
./diffusion_tiled_HT3_host
./diffusion_tiled_HT3_host
./diffusion_tiled_HT3_host

./diffusion_tiled_HT3_Large_host
./diffusion_tiled_HT3_Large_host
./diffusion_tiled_HT3_Large_host
else
echo "Please export MyNumberOfCores=(your number of cores here)"
fi


