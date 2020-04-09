#!/bin/bash

#script to run parameter convergence study

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

# run parameter convergence study with variety of cores
diffusion=(0.000001 0.1 0.2 0.3 0.4 0.5 0.6)
reaction=(0.000001 0.025 0.05 0.075 0.1 0.125 0.15)

for d in "${diffusion[@]}"; do
	for r in "${reaction[@]}"; do
		echo "submit diffusion: $d, reaction: $r"
		bsub -o "parameter-$d-$r" -n 96 -R fullnode -W 24:00 mpirun ./Performance-FisherSolver --diffCoef1 "$d" --reactCoef "$r" --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --name "parameter" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
	done
done
