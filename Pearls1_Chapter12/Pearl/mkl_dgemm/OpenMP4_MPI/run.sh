#!/bin/bash

runs=6 # for statistics
phi="7xxx"
phiNumCores=60 # for 7xxx for groupSize 16; set to 14 for 5xxx and 13 for 3xxx Xeon Phis
maxGroups=15   # for 7xxx for groupSize 16; set to 14 for 5xxx and 13 for 3xxx Xeon Phis
mpss="3.1"
compiler="intel14.0.3"
mpi="intelMPI4.1.1.036"
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
for pinning in "schedSetAffinity"
do
    cf_a=0.0
    cf_b=0.0
    cf_c=0.0
    cf=0.0
    echo "# groups, group size, performance [Gflops/s]" > mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${pinning}_${postfix}.log
    for g in 1 2 4 8 ${maxGroups}
    do
	# build program
	if [[ "$pinning" == "schedSetAffinity" ]]
	then
	    echo "# compile with schedSetAffinity"
	    mpiicpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION=${cf} -DPINNING -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    ompString="MIC_OMP_PLACES=\"\""
	    kmpString="MIC_KMP_AFFINITY=\"\""
	elif [[ "$pinning" == "compact" ]]
	then
	    echo "# compile with compact"
	    mpiicpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION=${cf} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    ompString="MIC_OMP_PLACES=\"\""
	    kmpString="MIC_KMP_AFFINITY=compact"
	elif [[ "$pinning" == "scatter" ]]
	then
	    echo "# compile with scatter"
	    mpiicpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION=${cf} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    ompString="MIC_OMP_PLACES=\"\""
	    kmpString="MIC_KMP_AFFINITY=scatter"
	else
	    echo "# compile"
	    mpiicpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION=${cf} -lrt -o mkl_dgemm.x mkl_dgemm.cpp
	    ompString="MIC_OMP_PLACES=\"\""
	    kmpString="MIC_KMP_AFFINITY=\"\""
	fi

	# execute and delete executable
	mpiString="-genv ${ompString} -genv ${kmpString}"
	for i in `seq 0 1 $(($g-1))`
	do
	    if [[ "$i" == "0" ]]
	    then
		mpiString="${mpiString} -np 1 -env MIC_KMP_PLACE_THREADS=$((${phiNumCores}/${g}))c,4t,0O ./mkl_dgemm.x"
	    else
		mpiString="${mpiString} : -np 1 -env MIC_KMP_PLACE_THREADS=$((${phiNumCores}/${g}))c,4t,$(($i*(${phiNumCores}/${g})))O ./mkl_dgemm.x"
	    fi
	done

	rm -f temp.log
	for r in `seq 1 1 $runs`
	do
	    sleep 2
	    mpirun ${mpiString} 2>&1 | tee temp_1.log
	    tail -n 1 temp_1.log >> temp.log
	done
	./average.x temp.log >> mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${pinning}_${postfix}.log
	rm -f mkl_dgemm.x temp_1.log temp.log

    done
done

# data transfer setup 1
cf_a=1.0
cf_b=0.0
cf_c=1.0
pinning="schedSetAffinity"
echo "# groups, group size, performance [Gflops/s]" > mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${pinning}_${postfix}.log
for g in 1 2 4 8 ${maxGroups}
do
    # build program
    echo "# compile with schedSetAffinity"
    mpiicpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION=${cf} -DCOPY_FRACTION_A=${cf_a} -DCOPY_FRACTION_B=${cf_b} -DCOPY_FRACTION_C=${cf_c} -DPINNING -lrt -o mkl_dgemm.x mkl_dgemm.cpp
    ompString="MIC_OMP_PLACES=\"\""
    kmpString="MIC_KMP_AFFINITY=\"\""

    # execute and delete executable
    mpiString="-genv ${ompString} -genv ${kmpString}"
    for i in `seq 0 1 $(($g-1))`
    do
	if [[ "$i" == "0" ]]
	then
	    mpiString="${mpiString} -np 1 -env MIC_KMP_PLACE_THREADS=$((${phiNumCores}/${g}))c,4t,0O ./mkl_dgemm.x"
	else
	    mpiString="${mpiString} : -np 1 -env MIC_KMP_PLACE_THREADS=$((${phiNumCores}/${g}))c,4t,$(($i*(${phiNumCores}/${g})))O ./mkl_dgemm.x"
	fi
    done
    
    rm -f temp.log
    for r in `seq 1 1 $runs`
    do
	sleep 2
	mpirun ${mpiString} 2>&1 | tee temp_1.log
	tail -n 1 temp_1.log >> temp.log
    done
    ./average.x temp.log >> mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${pinning}_${postfix}.log
    rm -f mkl_dgemm.x temp_1.log temp.log

done

# data transfer setup 2
cf_a=1.0
cf_b=1.0
cf_c=1.0
pinning="schedSetAffinity"
echo "# groups, group size, performance [Gflops/s]" > mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${pinning}_${postfix}.log
for g in 1 2 4 8 ${maxGroups}
do
    # build program
    echo "# compile with schedSetAffinity"
    mpiicpc -O3 -openmp -mkl -std=c++11 -D${output} -DMATRIX_SIZE=${size} -DGROUP_SIZE=${groupSize} -DCOPY_FRACTION=${cf} -DCOPY_FRACTION_A=${cf_a} -DCOPY_FRACTION_B=${cf_b} -DCOPY_FRACTION_C=${cf_c} -DPINNING -lrt -o mkl_dgemm.x mkl_dgemm.cpp
    ompString="MIC_OMP_PLACES=\"\""
    kmpString="MIC_KMP_AFFINITY=\"\""

    # execute and delete executable
    mpiString="-genv ${ompString} -genv ${kmpString}"
    for i in `seq 0 1 $(($g-1))`
    do
	if [[ "$i" == "0" ]]
	then
	    mpiString="${mpiString} -np 1 -env MIC_KMP_PLACE_THREADS=$((${phiNumCores}/${g}))c,4t,0O ./mkl_dgemm.x"
	else
	    mpiString="${mpiString} : -np 1 -env MIC_KMP_PLACE_THREADS=$((${phiNumCores}/${g}))c,4t,$(($i*(${phiNumCores}/${g})))O ./mkl_dgemm.x"
	fi
    done
    
    rm -f temp.log
    for r in `seq 1 1 $runs`
    do
	sleep 2
	mpirun ${mpiString} 2>&1 | tee temp_1.log
	tail -n 1 temp_1.log >> temp.log
    done
    ./average.x temp.log >> mkl_dgemm_performance_size${size}_groupSize${groupSize}_copyFraction${cf_a}_${cf_b}_${cf_c}_xeonPhi${phi}_mpss${mpss}_compiler${compiler}_${mpi}_${pinning}_${postfix}.log
    rm -f mkl_dgemm.x temp_1.log temp.log

done

rm -f average.x
