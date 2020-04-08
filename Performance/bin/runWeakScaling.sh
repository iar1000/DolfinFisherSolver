#!/bin/bash

#script to run weak scaling

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

# run weak scaling test on different mesh sizes
bsub -n 2 mpirun ./Performance-FisherSolver --name "weakscaling-40k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-80k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 16 mpirun ./Performance-FisherSolver --name "weakscaling-40k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 120 -R fullnode mpirun ./Performance-FisherSolver --name "weakscaling-40k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-4700k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"

bsub -n 1 mpirun ./Performance-FisherSolver --name "weakscaling-80k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-80k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 8 mpirun ./Performance-FisherSolver --name "weakscaling-80k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 60 mpirun ./Performance-FisherSolver --name "weakscaling-80k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-4700k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 463 mpirun ./Performance-FisherSolver --name "weakscaling-80k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-37000k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"

bsub -n 4 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 30 mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-4700k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 240 -R fullnode mpirun ./Performance-FisherSolver --name "weakscaling-150k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-37000k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"

bsub -n 1 mpirun ./Performance-FisherSolver --name "weakscaling-600k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-600k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 8 mpirun ./Performance-FisherSolver --name "weakscaling-600k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-4700k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
bsub -n 60 mpirun ./Performance-FisherSolver --name "weakscaling-600k" --meshname "lh-white-hull-flood-0-1-merge-5-dof-37000k.xml" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
