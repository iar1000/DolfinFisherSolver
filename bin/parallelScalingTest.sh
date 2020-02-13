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

echo "STRONG SCALING 24 PROCESSORS:"

PROCS=24

FILE_NAME="MPI12OMP2"
export OMP_NUM_THREADS=2
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				1 1 0.00001 1"

FILE_NAME="MPI6OMP4"
export OMP_NUM_THREADS=4
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				1 1 0.00001 1"

FILE_NAME="MPI4OMP6"				
export OMP_NUM_THREADS=6
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				1 1 0.00001 1"

FILE_NAME="MPI3OMP8"
export OMP_NUM_THREADS=8
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				1 1 0.00001 1"

FILE_NAME="MPI2OMP12"
export OMP_NUM_THREADS=12
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				1 1 0.00001 1"
				
