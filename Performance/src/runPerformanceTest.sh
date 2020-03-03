#!/bin/bash

#script to automatically run performance tests for 1-96 cores

SCRATCH=10000
MEM=12800

bsub mpirun ./PerformanceTest --ndofs 250000
PROCS=2
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=4
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=8
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=12
bsub -n "$PROCS" -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=24
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=36
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=48
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=72
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
PROCS=96
bsub -n "$PROCS" -R fullnode -R "rusage[scratch=$SCRATCH, mem=$(($MEM/$PROCS))]" mpirun ./PerformanceTest --ndofs 250000
