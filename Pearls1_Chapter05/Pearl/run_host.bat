: run_host.bat
if defined NUMBER_OF_CORES goto ncok
@ECHO NUMBER_OF_CORES not defined
@set /A NUMBER_OF_CORES=%NUMBER_OF_PROCESSORS% / 2
:ncok
@SET NUMBER_OF_CORES
@echo NUMBER_OF_THREADS==%NUMBER_OF_PROCESSORS%

set OMP_NUM_THREADS=
set KMP_AFFINITY=
set KMP_PLACE_THREADS=
diffusion_base_host
diffusion_base_host
diffusion_base_host

diffusion_base_Large_host
diffusion_base_Large_host
diffusion_base_Large_host

set KMP_AFFINITY=scatter
set OMP_NUM_THREADS=%NUMBER_OF_CORES%
diffusion_omp_host
diffusion_omp_host
diffusion_omp_host

diffusion_omp_Large_host
diffusion_omp_Large_host
diffusion_omp_Large_host

set KMP_AFFINITY=compact
set OMP_NUM_THREADS=
diffusion_omp_host
diffusion_omp_host
diffusion_omp_host

diffusion_omp_Large_host
diffusion_omp_Large_host
diffusion_omp_Large_host

set KMP_AFFINITY=scatter
set OMP_NUM_THREADS=%NUMBER_OF_CORES%
diffusion_ompvect_host
diffusion_ompvect_host
diffusion_ompvect_host

diffusion_ompvect_Large_host
diffusion_ompvect_Large_host
diffusion_ompvect_Large_host

set KMP_AFFINITY=compact
set OMP_NUM_THREADS=
diffusion_ompvect_host
diffusion_ompvect_host
diffusion_ompvect_host

diffusion_ompvect_Large_host
diffusion_ompvect_Large_host
diffusion_ompvect_Large_host

set KMP_AFFINITY=scatter
set OMP_NUM_THREADS=%NUMBER_OF_CORES%
diffusion_peel_host
diffusion_peel_host
diffusion_peel_host

diffusion_peel_Large_host
diffusion_peel_Large_host
diffusion_peel_Large_host

set KMP_AFFINITY=compact
set OMP_NUM_THREADS=
diffusion_peel_host
diffusion_peel_host
diffusion_peel_host

diffusion_peel_Large_host
diffusion_peel_Large_host
diffusion_peel_Large_host

set KMP_AFFINITY=scatter
set OMP_NUM_THREADS=%NUMBER_OF_CORES%
diffusion_tiled_host
diffusion_tiled_host
diffusion_tiled_host

diffusion_tiled_Large_host
diffusion_tiled_Large_host
diffusion_tiled_Large_host

set KMP_AFFINITY=compact
set OMP_NUM_THREADS=
diffusion_tiled_host
diffusion_tiled_host
diffusion_tiled_host

diffusion_tiled_Large_host
diffusion_tiled_Large_host
diffusion_tiled_Large_host

set KMP_AFFINITY=compact,fine
set OMP_NUM_THREADS=
diffusion_tiled_HT1_host
diffusion_tiled_HT1_host
diffusion_tiled_HT1_host

diffusion_tiled_HT1_Large_host
diffusion_tiled_HT1_Large_host
diffusion_tiled_HT1_Large_host

set KMP_AFFINITY=compact,fine
set OMP_NUM_THREADS=
diffusion_tiled_HT2_host
diffusion_tiled_HT2_host
diffusion_tiled_HT2_host

diffusion_tiled_HT2_Large_host
diffusion_tiled_HT2_Large_host
diffusion_tiled_HT2_Large_host

set KMP_AFFINITY=compact,fine
set OMP_NUM_THREADS=
diffusion_tiled_HT3_host
diffusion_tiled_HT3_host
diffusion_tiled_HT3_host

diffusion_tiled_HT3_Large_host
diffusion_tiled_HT3_Large_host
diffusion_tiled_HT3_Large_host

