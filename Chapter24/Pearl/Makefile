SIZE=4000
TRIALS=32
CPU_CORES=`cat /proc/cpuinfo | grep "processor" | wc -l`
MIC_CORES=`micinfo -d 0 | grep "Active Cores" | awk '{print $$7}'`

all: 
	cd 1; make
	cd 2; make
	cd 3; make
	cd 4; make
	cd 5; make

run-cpu:
	@cd 1; OMP_NUM_THREADS=$(CPU_CORES) KMP_AFFINITY=compact ./runme-CPU $(SIZE) $(TRIALS); echo
	@cd 2; OMP_NUM_THREADS=$(CPU_CORES) KMP_AFFINITY=compact ./runme-CPU $(SIZE) $(TRIALS); echo
	@cd 3; OMP_NUM_THREADS=$(CPU_CORES) KMP_AFFINITY=compact ./runme-CPU $(SIZE) $(TRIALS); echo
	@cd 4; OMP_NUM_THREADS=$(CPU_CORES) KMP_AFFINITY=compact ./runme-CPU $(SIZE) $(TRIALS); echo
	@cd 5; OMP_NUM_THREADS=$(CPU_CORES) KMP_AFFINITY=compact ./runme-CPU $(SIZE) $(TRIALS); echo

run-mic:
	@cd 1; scp -q ./runme-MIC mic0:~/; ssh mic0 KMP_PLACE_THREADS=$(MIC_CORES)C,4t ~/runme-MIC $(SIZE) $(TRIALS); echo
	@cd 2; scp -q ./runme-MIC mic0:~/; ssh mic0 KMP_PLACE_THREADS=$(MIC_CORES)C,4t ~/runme-MIC $(SIZE) $(TRIALS); echo
	@cd 3; scp -q ./runme-MIC mic0:~/; ssh mic0 KMP_PLACE_THREADS=$(MIC_CORES)C,4t ~/runme-MIC $(SIZE) $(TRIALS); echo
	@cd 4; scp -q ./runme-MIC mic0:~/; ssh mic0 KMP_PLACE_THREADS=$(MIC_CORES)C,4t ~/runme-MIC $(SIZE) $(TRIALS); echo
	@cd 5; scp -q ./runme-MIC mic0:~/; ssh mic0 KMP_PLACE_THREADS=$(MIC_CORES)C,4t ~/runme-MIC $(SIZE) $(TRIALS); echo

clean:
	cd 1; make clean
	cd 2; make clean
	cd 3; make clean
	cd 4; make clean
	cd 5; make clean
