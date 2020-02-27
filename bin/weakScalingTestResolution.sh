#!/bin/bash

# bash script to run weak scaling test for FisherSolver
# weak scaling is used to measure distributed memory scaling (MPI)
# this test suite is developed for euler architecture, 24 cores (48 virtual)


SCRATCH=10000
MEM=12800
TIMESTEPS=100000

echo "RUNNING WEAK SCALING TEST (RESOLUTION) FOR FISHER-SOLVER"
echo "   SCRATCH=$SCRATCH MB"
echo "   MEMORY=$MEM MB"
				
## MESH RESOLUTION
# BASELINE 
PROCS=1
MESH_NAME=rect-10on10-res-25.h5
FILE_NAME="BASELINE-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]"  ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=2
MESH_NAME=rect-10on10-res-35.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
PROCS=4
MESH_NAME=rect-10on10-res-50.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=8
MESH_NAME=rect-10on10-res-70.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=12
MESH_NAME=rect-10on10-res-85.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

PROCS=24
MESH_NAME=rect-10on10-res-120.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
				
PROCS=36
MESH_NAME=rect-10on10-res-140.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1
					
PROCS=48
MESH_NAME=rect-10on10-res-170.h5
FILE_NAME="${PROCS}-$PP{RESOLUTION}WEAK"
bsub -n "$PROCS" -o "${PROCS}WEAK-out.txt" -W 08:00 -R "rusage[scratch=10000, mem=$(($MEM/$PROCS))]" -R fullnode mpirun -n "$PROCS" ./FisherSolver \
				"../output" "../mesh" "$MESH_NAME" "WEAKSCALING" "$FILE_NAME" "${PROCS}iterationdata"\
				1 1 1 1 1 \
				0.000000001 0.000001 1 0.1 0 \
				0.013 0.0013 0.025 1\
				2 1 0.00001 1

				
