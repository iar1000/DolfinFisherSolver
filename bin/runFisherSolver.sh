#!/bin/bash

# bash script to run the FisherSolver

# Get Simulation configuration data
# read there for further information
source FisherSolver.config
# other parameters default values
VERBOSE=2
FOLDER_NAME=$(date +%F)
FILE_NAME="out"

# default values can be overwritten by following command line arguments
# --verbose: set verbosity level
#		1: no simulation output
#		2: only simulation progress (default)
#		3: progress and information about each iteration
# --outputfolder: name of the output subfolder, must be valid directory name otherwise undefined behavior
# 		default: current date
# --outputfile: name of output files
#		of output is always .pvd
#		default: "out"
# --diffusionwhite: diffusion parameter of white matter
# --diffusiongrey: diffusion parameter of grey matter
# --rho: reaction parameter
while [[ "$#" -gt 0 ]]; do case $1 in
  -v|--verbose) VERBOSE="$2"; shift;;
  -ofo|--outputfolder) FOLDER_NAME="$2"; shift;;
  -ofi|--outputfile) FILE_NAME="$2"; shift;;
  -dw|--diffusionwhite) DIFFUSION_W="$2"; shift;;
  -dg|--diffusiongrey) DIFFUSION_G="$2"; shift;;
  -rh|--rho) RHO="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done


echo "RUNNING FISHER-SOLVER"
echo "PARAMETERS:"
echo "  VERBOSE= $VERBOSE"
echo "  OUTPUT FOLDER= output/$FOLDER_NAME"
echo "  OUTPUT FILE= $FILE_NAME"
echo "  MESH NAME= $MESH_NAME"
echo "  DIMENSIONALITY= $SPATIAL_DIMENSION"
echo "  INITIALIZATION @ ($CX, $CY, $CZ) VALUE $VALUE"
echo "  DIFFUSION PARAMETERS= $DIFFUSION_W (W), $DIFFUSION_G (G)"
echo "  REACTION PARAMETERS = $RHO"


