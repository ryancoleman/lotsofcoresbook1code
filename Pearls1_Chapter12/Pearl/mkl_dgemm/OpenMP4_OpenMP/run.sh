#!/bin/bash

runs=6 # for statistics
phi="7xxx"
maxGroups=15  # for 7xxx; set to 14 for 5xxx and 13 for 3xxx Xeon Phis
mpss="3.1"
compiler="intel14.0.3"
postfix="XXX"
size="2048L"
groupSize="16"
#output="BATCH_OUTPUT"
output=""

icpc -o average.x average.cpp

# multi-threaded MKL: will be modified in program
export MKL_NUM_THREADS=2
export MIC_MKL_NUM_THREADS=2

# mic env-prefix
export MIC_ENV_PREFIX=MIC

# performance for different pinning schemes (no data transfers)
for pinning in "schedSetAffinity" "fixedCompact" "fixedScatter" "fixedCores" "compact" "scatter" "cores" "noAffinity"
do
    cf_a=0.0
    cf_b=0.0
    cf_c=0.0
    cf=0.0
    echo "# groups, group size, performance [Gflops/s]" > mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${pinning}_${postfix}.log
    for g in 1 2 4 8 ${maxGroups}
    do
	# build program
	if [[ "$pinning" == "schedSetAffinity" ]]
	then
	    echo "# compile with schedSetAffinity"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -DPINNING -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	elif [[ "$pinning" == "fixedCompact" ]]
	then
	    echo "# compile with fixedCompact"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -DFIX -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=
	    export MIC_KMP_AFFINITY=compact
	elif [[ "$pinning" == "fixedScatter" ]]
	then
	    echo "# compile with fixedScatter"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -DFIX -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=
	    export MIC_KMP_AFFINITY=scatter
	elif [[ "$pinning" == "fixedCores" ]]
	then
	    echo "# compile with fixedCores"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -DFIX -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=cores
	    export MIC_KMP_AFFINITY=
	elif [[ "$pinning" == "compact" ]]
	then
	    echo "# compile with compact"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=
	    export MIC_KMP_AFFINITY=compact
	elif [[ "$pinning" == "scatter" ]]
	then
	    echo "# compile with scatter"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=
	    export MIC_KMP_AFFINITY=scatter
	elif [[ "$pinning" == "cores" ]]
	then
	    echo "# compile with cores"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=cores
	    export MIC_KMP_AFFINITY=
	else
	    echo "# compile"
	    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    export MIC_OMP_PLACES=
	    export MIC_KMP_AFFINITY=
	fi

	# execute and delete executable
	rm -f temp.log
	for r in `seq 1 1 $runs`
	do
	    sleep 5
	    ./mkl_dgemm.x 2>&1 | tee temp_1.log
	    tail -n 1 temp_1.log >> temp.log
	done
	./average.x temp.log >> mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${pinning}_${postfix}.log
	rm -f mkl_dgemm.x temp_1.log temp.log
    done
done

# data transfer setup 1
cf_a=1.0
cf_b=0.0
cf_c=1.0
pinning="schedSetAffinity"
echo "# groups, group size, performance [Gflops/s]" > mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${pinning}_${postfix}.log
for g in 1 2 4 8 ${maxGroups}
do
    # build program
    echo "# compile with schedSetAffinity"
    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION_A=${cf_a} -DCOPY_FRACTION_B=${cf_b} -DCOPY_FRACTION_C=${cf_c} -DPINNING -lrt -o mkl_dgemm.x mkl_dgemm.cpp

    # execute and delete executable
    rm -f temp.log
    for r in `seq 1 1 $runs`
    do
	sleep 5
	./mkl_dgemm.x 2>&1 | tee temp_1.log
	tail -n 1 temp_1.log >> temp.log
    done
    ./average.x temp.log >> mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${pinning}_${postfix}.log
    rm -f mkl_dgemm.x temp_1.log temp.log
done

# data transfer setup 2
cf_a=1.0
cf_b=1.0
cf_c=1.0
pinning="schedSetAffinity"
echo "# groups, group size, performance [Gflops/s]" > mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${pinning}_${postfix}.log
for g in 1 2 4 8 ${maxGroups}
do
    # build program
    echo "# compile with schedSetAffinity"
    icpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DNUM_GROUPS=${g} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION_A=${cf_a} -DCOPY_FRACTION_B=${cf_b} -DCOPY_FRACTION_C=${cf_c} -DPINNING -lrt -o mkl_dgemm.x mkl_dgemm.cpp

    # execute and delete executable
    rm -f temp.log
    for r in `seq 1 1 $runs`
    do
	sleep 5
	./mkl_dgemm.x 2>&1 | tee temp_1.log
	tail -n 1 temp_1.log >> temp.log
    done
    ./average.x temp.log >> mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${pinning}_${postfix}.log
    rm -f mkl_dgemm.x temp_1.log temp.log
done

rm -f average.x