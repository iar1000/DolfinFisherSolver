#!/bin/bash

#script to automatically run performance tests for 1-96 cores

SCRATCH=10000
MEM=750000
NDOFS=500000
TOL=0.00000001
TYPE=1

# overwrite default by command line arguments
while [[ "$#" -gt 0 ]]; do case $1 in
  -dof|--ndofs) NDOFS="$2"; shift;;
  -me|--memory) NDOFS="$2"; shift;;
  -tol|--tolerance) TOL="$2"; shift;;
  -t|--type) TYPE="$2"; shift;;
 
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

bsub -R "rusage[scratch=$SCRATCH, mem=50000]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=2
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=4
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=8
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=12
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=24
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=36
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=48
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=72
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
PROCS=96
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE"
