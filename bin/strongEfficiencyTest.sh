#!/bin/bash

# bash script to run the scaling tests for FisherSolver
# this test suite is developed for euler architecture, 24 cores (48 virtual)

# Get Mesh of desired production run
# rest of the parameters are hardcoded default values
source FisherSolver.config

SCRATCH=10000
MEM=12800
TIMESTEPS=10000 # @IMPORANT: Since bash doesn't include floating point arithmetic, must set dt yourself or change $TIME
TIME=1			# per T => 10'000 steps of size 0.0001

echo "RUNNING PARALLEL EFFICIENCY TESTS FOR FISHER-SOLVER"
echo "   SCRATCH=$SCRATCH MB"
echo "   MEMORY=$MEM MB"

##########################################################
# BASELINE #############################
#########################################################
PROCS=1
echo "STRONG SCALING $PROCS PROCESSORS:"

FILE_NAME="BASELINE-${TIMESTEPS}STRONG"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "PARALLELEFFICIENCYTEST" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$TIME" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

##########################################################
# 4 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=4
echo "STRONG SCALING $PROCS PROCESSORS:"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "mpirun -n $PROCS --report-bindings ./FisherSolver \
				../output ../mesh $MESH_NAME PARALLELEFFICIENCYTEST $FILE_NAME \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"
				

##########################################################
# 8 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=8
echo "STRONG SCALING $PROCS PROCESSORS:"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "mpirun -n $PROCS --report-bindings ./FisherSolver \
				../output ../mesh $MESH_NAME PARALLELEFFICIENCYTEST $FILE_NAME \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"
				

##########################################################
# 12 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=12
echo "STRONG SCALING $PROCS PROCESSORS:"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "mpirun -n $PROCS --report-bindings ./FisherSolver \
				../output ../mesh $MESH_NAME PARALLELEFFICIENCYTEST $FILE_NAME \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"

##########################################################
# 24 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=24
echo "STRONG SCALING $PROCS PROCESSORS:"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG-MPI24OMP1"
export OMP_NUM_THREADS=1
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS --report-bindings ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG-MPI12OMP2"
export OMP_NUM_THREADS=2
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS --report-bindings ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG-MPI6OMP4"
export OMP_NUM_THREADS=4
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS --report-bindings ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG-MPI4OMP6"				
export OMP_NUM_THREADS=6
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS --report-bindings ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG-MPI3OMP8"
export OMP_NUM_THREADS=8
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS --report-bindings ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"

FILE_NAME="${PROCS}-${TIMESTEPS}STRONG-MPI2OMP12"
export OMP_NUM_THREADS=12
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" "unset LSB_AFFINITY_HOSTFILE ; mpirun -n $(($PROCS/$OMP_NUM_THREADS)) --map-by node:PE=$OMP_NUM_THREADS --report-bindings ./FisherSolver \
				\"../output\" \"../mesh\" \"$MESH_NAME\" \"PARALLELEFFICIENCYTEST\" \"$FILE_NAME\" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 $TIME 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1"
				
