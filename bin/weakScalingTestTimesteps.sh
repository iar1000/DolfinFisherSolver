#!/bin/bash

# bash script to run weak scaling test for FisherSolver
# weak scaling is used to measure distributed memory scaling (MPI)
# this test suite is developed for euler architecture, 24 cores (48 virtual)


SCRATCH=10000
MEM=12800
TIMESTEPS=5000 # per processor

echo "RUNNING WEAK SCALING TEST (TIMESTEPS) FOR FISHER-SOLVER"
echo "   SCRATCH=$SCRATCH MB"
echo "   MEMORY=$MEM MB"

## TIMESTEPS
MESH_NAME=rect-10on10-res-100.h5
# BASELINE 
PROCS=1
FILE_NAME="BASELINE-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"  ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.005 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=2
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.01 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
PROCS=4
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.02 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=8
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.04 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=12
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.06 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=24
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.12 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
PROCS=36
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.18 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
					
PROCS=48
FILE_NAME="${PROCS}-$PP{TIMESTEPS}WEAK"
bsub -n "$PROCS" -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.24 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1