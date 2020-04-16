#!/bin/bash

# bash script to run the FisherSolver

# Get Simulation configuration data
# read there for further information
source FisherSolver.config

# default values and argument description
# --mpiprocess: set number of mpi processes
MPI_PROCESS=2
# --verbose: set verbosity level
#		1: no simulation output
#		2: only simulation progress (default)
#		3: progress and information about each iteration
#		4: debug
VERBOSE=4 
# --outputfolder: name of the output subfolder, must be valid directory name otherwise undefined behavior
FOLDER_NAME=$(date +%F)
# --outputfile: name of output files
#		datatype of output is always .pvd
FILE_NAME="out"
# --csvfile: name of csv file for iteration data
#		datatype of iteration data is always .csv
CSV_NAME="iterationdata"
# --framerate: number of frames saved to file per 1 time unit (TEND=10, FRAMERATE=30 => 300 FRAMES SAVED)
FRAMERATE=30

OUTPUT_PARENTFOLDER="../output" 		# path to output parent folder
MESH_PARENTFOLDER="../mesh" 	# path to mesh parent folder

# overwrite default by command line arguments
while [[ "$#" -gt 0 ]]; do case $1 in
  -n|--mpiprocess) MPI_PROCESS="$2"; shift;;
  -v|--verbose) VERBOSE="$2"; shift;;
  -ofo|--outputfolder) FOLDER_NAME="$2"; shift;;
  -ofi|--outputfile) FILE_NAME="$2"; shift;;
  -csv|--csvfile) CSV_NAME="$2"; shift;;
  -fr|--framerate) FRAMERATE="$2"; shift;;
  # overwrite .config
  -m|--mesh) MESH_NAME="$2"; shift;;
  -dw|--diffusionwhite) DIFFUSION_W="$2"; shift;;
  -dg|--diffusiongrey) DIFFUSION_G="$2"; shift;;
  -rh|--rho) RHO="$2"; shift;;
  -th|--theta) THETA="$2"; shift;;
  -dt|--dtstart) DTSTART="$2"; shift;;
  -T|--Tend) TEND="$2"; shift;;
  -type|--dttype) TIMEADAPTION="$2"; shift;;
  -rs|--richsafe) RICHARDSONSAFETY="$2"; shift;;
  -rt|--richtol) RICHARDSONTOL="$2"; shift;;
  -rl|--runlength) RUNLENGTH="$2"; shift;;
  -ls|--solver) KRYLOVSOLVER="$2"; shift;;
  -pc|--preconditioner) KRYLOVPREC="$2"; shift;;
    
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done


echo "RUNNING FISHER-SOLVER WITH $MPI_PROCESS PROCESS"
echo "PARAMETERS:"
echo "  VERBOSE= $VERBOSE"
echo "  OUTPUT PATH= $OUTPUT_PARENTFOLDER/$FOLDER_NAME"
echo "  OUTPUT FILE= $FILE_NAME (FRAME RATE= $FRAMERATE)"
echo "  CSV FILE= $CSV_NAME"
echo "  MESH PATH= $MESH_PARENTFOLDER/$MESH_NAME"
echo "	CALCULATE TRANSLATION AND STOP= $CALCTRANS"
echo "  INITIALIZATION @ ($CX, $CY, $CZ) VALUE $VALUE"
echo "  DIFFUSION PARAMETERS= $DIFFUSION_W (W), $DIFFUSION_G (G)"
echo "  REACTION PARAMETERS = $RHO"
echo "  DISCRETIZATION PARAMETER = $THETA"
echo "  TIMESTEPPING: to T=$TEND with dt=$DTSTART"
echo "  TIMEADAPTION: $TIMEADAPTION (tol= $RICHARDSONTOL, safety= $RICHARDSONSAFETY)"
echo "  RUNLENGTH: $RUNLENGTH"
echo "  LINEAR SOLVER: $KRYLOVSOLVER, PRECONDITIONER: $KRYLOVPREC"

mpirun -n "$MPI_PROCESS" ./FisherSolver \
				"$OUTPUT_PARENTFOLDER" 	 "$MESH_PARENTFOLDER" "$MESH_NAME" "$FOLDER_NAME" "$FILE_NAME" "$CSV_NAME" \
				"$CX" "$CY" "$CZ" "$RADIUS" "$VALUE" \
				"$DTMIN" "$DTSTART" "$DTMAX" "$TEND" "$FRAMERATE" \
				"$DIFFUSION_W" "$DIFFUSION_G" "$RHO" "$THETA"\
				"$VERBOSE" "$TIMEADAPTION" "$RICHARDSONTOL" "$RICHARDSONSAFETY" "$RUNLENGTH"\
				"$NEWTONRESIDUALTOLREL" "$NEWTONRESIDUALTOLABS" "$KRYLOVRESIDUALTOLREL" "$KRYLOVRESIDUALTOLABS" \
				"$KRYLOVSOLVER" "$KRYLOVPREC" \
				"$CALCTRANS" "$TRANSX" "$TRANSY" "$TRANSZ"
				
				
