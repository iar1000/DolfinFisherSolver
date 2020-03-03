#!/bin/bash

#script to automatically run performance tests for 1-96 cores

SCRATCH=10000
MEM=3000000
NDOFS=500000

bsub mpirun ./PerformanceTest --ndofs 250000
PROCS=2
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=4
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=8
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=12
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=24
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=36
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=48
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=72
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
PROCS=96
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs "$NDOFS"
