#!/bin/bash

# bash script to run weak scaling test for FisherSolver
# weak scaling is used to measure distributed memory scaling (MPI)
# this test suite is developed for euler architecture, 24 cores (48 virtual)


MESH_NAME=rect-100on100-res-100.h5

SCRATCH=10000
MEM=12800
TIME=0.5			# per T => 10'000 steps of size 0.0001

echo "RUNNING WEAK TEST FOR FISHER-SOLVER"
echo "   SCRATCH=$SCRATCH MB"
echo "   MEMORY=$MEM MB"

##########################################################
# BASELINE #############################
#########################################################
PROCS=1
FILE_NAME="BASELINE-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$TIME" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

##########################################################
# 2 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=2
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
##########################################################
# 4 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=4
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

##########################################################
# 8 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=8
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1


##########################################################
# 12 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=12
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

##########################################################
# 24 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=24
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
##########################################################
# 36 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=36
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
##########################################################
# 48 PROCESSES SCALING TESTS #############################
#########################################################
PROCS=48
FILE_NAME="${PROCS}-${TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"-R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" \
				1 1 1 1 1 \
				0.000000001 0.0001 1 "$(($TIME*$PROCS))" 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

				
