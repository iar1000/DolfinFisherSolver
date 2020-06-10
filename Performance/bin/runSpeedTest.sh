#!/bin/bash

#script to run speed test with plial meshes mid- and big-sized

# type 1: all combinations of ls and pc
TYPE=3
TOL=0.00000001
KRYLNONZERO=0
MESHSIZE=5500
# overwrite default by command line arguments
while [[ "$#" -gt 0 ]]; do case $1 in
  -tol|--tolerance) TOL="$2"; shift;;
  -knz| --nonzero) KRYLNONZERO="$2"; shift;;
  -t| --type) TYPE="$2"; shift;;
  -m| --size) MESHSIZE="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# run speed test with variety of cores for different mesh sizes
if [ "$MESHSIZE" -eq 5500 ]
then
	fncores=(24 36 48 96 120 180 240 360)
elif [ "$MESHSIZE" -eq 40000 ]
then
	fncores=(120 180 240 360 480 540 600)
fi

echo "${cores[*]}"
echo "${fncores[*]}"

for c in "${cores[@]}"; do
	echo "submit $c core job"
	bsub -o "speed-${MESHSIZE}-$c" -n "$c" -W 08:00  mpirun ./Performance-FisherSolver --meshname "lh-plial-dof-${MESHSIZE}k.h5" --name "speed${MESHSIZE}" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
done

for c in "${fncores[@]}"; do
	echo "submit $c core job fullnode"
	bsub -o "speed-${MESHSIZE}-$c" -n "$c" -W 08:00 -R fullnode mpirun ./Performance-FisherSolver --meshname "lh-plial-dof-${MESHSIZE}k.h5" --name "speed${MESHSIZE}" --type "$TYPE" --newton_tol "$TOL" --krylovnonzero "$KRYLNONZERO"
done