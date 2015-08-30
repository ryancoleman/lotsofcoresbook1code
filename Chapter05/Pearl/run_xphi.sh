#!/bin/sh
unset OMP_NUM_THREADS
unset KMP_AFFINITY
unset OMP_PLACE_THREADS
./diffusion_base_xphi
./diffusion_base_xphi
./diffusion_base_xphi

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=180
./diffusion_omp_xphi
./diffusion_omp_xphi
./diffusion_omp_xphi

./diffusion_omp_Large_xphi
./diffusion_omp_Large_xphi
./diffusion_omp_Large_xphi

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=180
./diffusion_ompvect_xphi
./diffusion_ompvect_xphi
./diffusion_ompvect_xphi

./diffusion_ompvect_Large_xphi
./diffusion_ompvect_Large_xphi
./diffusion_ompvect_Large_xphi

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=180
./diffusion_peel_xphi
./diffusion_peel_xphi
./diffusion_peel_xphi

./diffusion_peel_Large_xphi
./diffusion_peel_Large_xphi
./diffusion_peel_Large_xphi

export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=180
./diffusion_tiled_xphi
./diffusion_tiled_xphi
./diffusion_tiled_xphi

./diffusion_tiled_Large_xphi
./diffusion_tiled_Large_xphi
./diffusion_tiled_Large_xphi

export KMP_AFFINITY=compact
unset OMP_NUM_THREADS
# (run with 240 threads)
./diffusion_tiled_HT1_xphi
./diffusion_tiled_HT1_xphi
./diffusion_tiled_HT1_xphi

./diffusion_tiled_HT1_Large_xphi
./diffusion_tiled_HT1_Large_xphi
./diffusion_tiled_HT1_Large_xphi

./diffusion_tiled_HT2_xphi
./diffusion_tiled_HT2_xphi
./diffusion_tiled_HT2_xphi

./diffusion_tiled_HT2_Large_xphi
./diffusion_tiled_HT2_Large_xphi
./diffusion_tiled_HT2_Large_xphi

./diffusion_tiled_HT3_xphi
./diffusion_tiled_HT3_xphi
./diffusion_tiled_HT3_xphi

./diffusion_tiled_HT3_Large_xphi
./diffusion_tiled_HT3_Large_xphi
./diffusion_tiled_HT3_Large_xphi


