#!/bin/bash

#script to run weak scaling for lh-plial meshes with 150k elements per core

# type 1: all combinations of ls and pc
# type 2: cg + hypre euclid
# type 4: bicgstab + hypre euclid
TYPE=2
TOL=0.00000001
KRYLNONZERO=0
# overwrite default by command line arguments
while [[ "$#" -gt 0 ]]; do case $1 in
  -tol|--tolerance) TOL="$2"; shift;;
  -knz| --nonzero) KRYLNONZERO="$2"; shift;;
  -t| --type) TYPE="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# 150k elements per core
bsub -o "weak-150-1" -n 1 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-plial-015mio.xdmf" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -o "weak-150-2" -n 2 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-plial-03mio.xdmf" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -o "weak-150-6" -n 6 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-plial-09mio.xdmf" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -o "weak-150-20" -n 20 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-plial-3mio.xdmf" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -o "weak-150-40" -n 40 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-plial-6mio.xdmf" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"

