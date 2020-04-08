#!/bin/bash

#script to run weak scaling of test

bsub -n 96 -R "rusage[scratch=10000, mem=10000]" mpirun ./PerformanceTest --ndofs 500000 --tol 0.00000001 --type 1 --krylovnonzero 1

ALL=1
TOL=0.00000001
KRYLNONZERO=0
# overwrite default by command line arguments
while [[ "$#" -gt 0 ]]; do case $1 in
  -tol|--tolerance) TOL="$2"; shift;;
  -knz| --nonzero) KRYLNONZERO="$2"; shift;;
  -a| --all) ALL="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

bsub -R "rusage[scratch=$SCRATCH, mem=50000]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=2
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=4
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=8
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=12
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=24
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=36
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=48
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=72
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
PROCS=96
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS" --tol "$TOL" --type "$TYPE" --krylovnonzero "$KRYLNONZERO"
