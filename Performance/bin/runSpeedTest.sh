#!/bin/bash

#script to run speed test

# type 1: all combinations of ls and pc
TYPE=1
TOL=0.00000001
KRYLNONZERO=0
# overwrite default by command line arguments
while [[ "$#" -gt 0 ]]; do case $1 in
  -tol|--tolerance) TOL="$2"; shift;;
  -knz| --nonzero) KRYLNONZERO="$2"; shift;;
  -t| --type) TYPE="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# run speed test with variety of cores
cores=(1 4 16)
fncores=(24 36 48 96 120 168 240 360)

for c in "${cores[@]}"; do
	echo "submit $c core job"
	bsub -n "$c" mpirun ./Performance-FisherSolver --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --name "speed" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
done

for c in "${fncores[@]}"; do
	echo "submit $c core job fullnode"
	bsub -n "$c" -R fullnode mpirun ./Performance-FisherSolver --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --name "speed" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
done