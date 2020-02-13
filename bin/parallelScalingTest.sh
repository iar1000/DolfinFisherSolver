#!/bin/bash

# bash script to run the scaling tests for FisherSolver
# this test suite is developed for euler architecture, 24 cores (48 virtual)

# Get Simulation configuration data
# read there for further information
# the same config file is used such that the scalability adapts to the desired parameters
source FisherSolver.config

SCRATCH=10000
MEM=12800

echo "RUNNING PARALLEL EFFICIENCY TESTS FOR FISHER-SOLVER"
echo "   SCRATCH=$SCRATCH MB"
echo "   MEMORY=$MEM MB"
echo "STRONG SCALING:"
echo "$OMP_NUM_THREADS"

export OMP_NUM_THREADS=4
bsub -I -n 24 -R fullnode -R "rusage[scratch=10000, mem=5000]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n 8 --map-by-node:PE=4 ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				1 1 0.00001 1"
				
