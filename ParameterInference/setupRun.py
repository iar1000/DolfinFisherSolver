# script to generate a run enviroment for parametric inference
# automatically creates a run directory in /runs with generated cases

import sys
import os
from datetime import date
import numpy as np
import stat


# https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros
def floatToString(inputValue):
    return ('%.15f' % inputValue).rstrip('0').rstrip('.')


# Parameter inference
diffusion_fac = 10  # factor D_w / D_g
diffusion_min = 0.1  # min D_W
diffusion_max = 0.6  # max D_W
diffusion_steps = 2  # discretization steps of parameter range
rho_min = 0.025  # min rho
rho_max = 0.25  # max rho
rho_steps = 2  # discretization steps of parameter range
radius_min = 1  # min r0
radius_max = 2  # max r0
radius_steps = 2  # discretization steps of parameter range
# Fixed Arguments
# Runtime arguments
mpiprocs = 2
verbosity = 3
# Model arguments
date = date.today().strftime("%d-%m-%Y")
initial_condition = [49, 130, 40, 3]  # x, y, z, radius
translation = [25, 21, 30]
# Solver arguments
mesh_name = "lh-white-hull-flood-0-1-merge-5-dof-80k.h5"
T_end = 0.001
dt = 0.0001
solver = "cg"
preconditioner = "hypre_euclid"
rich_tol = 0.000001
rich_safe = 0.95
newton_abs = 0.00000001
newton_rel = 0.00000001
framerate = -1

print('Generate run-case')
print("\tArguments: "
      "\n\t\t diffusion factor : {}"
      "\n\t\t D_w range : [{}, {}]"
      "\n\t\t rho range : [{}, {}]"
      "\n\t\t radius range : [{}, {}]"
      "\n\t\t discretization steps : "
      "\n\t\t\t D : {}"
      "\n\t\t\t rho : {}"
      "\n\t\t\t radius : {}".format(diffusion_fac, diffusion_min, diffusion_max, rho_min, rho_max, radius_min,
                                   radius_max, diffusion_steps, rho_steps, radius_steps))

# create parameter space
diffusion_range = diffusion_max - diffusion_min
diffusion_stepsize = diffusion_range / diffusion_steps
diffusion_pspace = np.arange(diffusion_min, diffusion_max + 0.01, diffusion_stepsize).tolist()
diffusion_pspace = ['%.3f' % elem for elem in diffusion_pspace]
rho_range = rho_max - rho_min
rho_stepsize = rho_range / rho_steps
rho_pspace = np.arange(rho_min, rho_max + 0.01, rho_stepsize).tolist()
rho_pspace = ['%.3f' % elem for elem in rho_pspace]
radius_range = radius_max - radius_min
radius_stepsize = radius_range / radius_steps
radius_pspace = np.arange(radius_min, radius_max + 0.01, radius_stepsize).tolist()
radius_pspace = ['%.3f' % elem for elem in radius_pspace]

print("Parameter spaces : \n"
      "\t D_w : {} \n"
      "\t rho : {}\n"
      "\t radius : {}\n"
      "\t Total cases : {}\n".format(diffusion_pspace, rho_pspace, radius_pspace,
                                     len(radius_pspace) * len(diffusion_pspace) * len(rho_pspace)))

# create run directory
rundir_path = "rundir-{}-minmax-D-{}_{}-rho-{}_{}-radius-{}_{}" \
    .format(date, diffusion_min, diffusion_max, rho_min, rho_max, radius_min, radius_max).replace('.', '')
parent_path = "runs/" + rundir_path
try:
    duplicate_counter = 0
    worker_path = parent_path
    while os.path.isdir(worker_path):
        duplicate_counter = duplicate_counter + 1
        worker_path = parent_path + "-{}".format(duplicate_counter)
    os.mkdir(worker_path)
    parent_path = worker_path
    print("created run directory {}".format(worker_path))
except Exception as e:
    print(e)
    quit()
print("run directory path {}".format(rundir_path))
print("parent path {}".format(parent_path))

# create case directories inside run directory
case_dirs = []  # relative path to case folder, d, r
for d in diffusion_pspace:
    for rho in rho_pspace:
        for radius in radius_pspace:
            case_path = "case-D-{}-rho-{}-radius-{}".format(d, rho, radius).replace(".", "")
            try:
                os.mkdir(parent_path + "/" + case_path)
                case_dirs.append([case_path, float(d), float(rho), float(radius)])
            except Exception as e:
                print("failed creating case directory {}".format(parent_path + case_path))
                print(e)
                quit()
print("created {} case directories".format(len(case_dirs)))

# update parameter spaces to be floats
diffusion_pspace = [float(i) for i in diffusion_pspace]
rho_pspace = [float(i) for i in rho_pspace]
radius_pspace = [float(i) for i in radius_pspace]

# create submission file for all cases
with open(parent_path + '/submit-all.sh', 'w') as submit_bash:
    command = '''\
#! /bin/bash
# automatically generated submission file
# the defined standart parameters can be overwritten by command line parameters

MPI_PROCESS={}
VERBOSE={}
FRAMERATE={}
MESH_NAME={}
DTSTART={}
TEND={}
TIMEADAPTION=2
RICHARDSONSAFETY={}
RICHARDSONTOL={}
KRYLOVSOLVER={}
KRYLOVPREC={}
NEWTONABS={}
NEWTONREL={}

while [[ "$#" -gt 0 ]]; do case $1 in
  -n|--mpiprocess) MPI_PROCESS="$2"; shift;;
  -v|--verbose) VERBOSE="$2"; shift;;
  -fr|--framerate) FRAMERATE="$2"; shift;;
  # overwrite .config
  -m|--mesh) MESH_NAME="$2"; shift;;
  -dt|--dtstart) DTSTART="$2"; shift;;
  -T|--Tend) TEND="$2"; shift;;
  -type|--dttype) TIMEADAPTION="$2"; shift;;
  -rs|--richsafe) RICHARDSONSAFETY="$2"; shift;;
  -rt|--richtol) RICHARDSONTOL="$2"; shift;;
  -ls|--solver) KRYLOVSOLVER="$2"; shift;;
  -pc|--preconditioner) KRYLOVPREC="$2"; shift;;  
  -nr|--newtonrel) NEWTONREL="$2"; shift;;  
  -na|--newtonabs) NEWTONABS="$2"; shift;;  
  
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# automatically added submissions
'''.format(mpiprocs, verbosity, framerate, mesh_name, dt, T_end, floatToString(rich_safe), floatToString(rich_tol),
           solver, preconditioner, floatToString(newton_abs),
           floatToString(newton_rel))
    submit_bash.write(command)
os.chmod(parent_path + '/submit-all.sh', 0o755)

# create submission file for corner cases
with open(parent_path + '/submit-corner.sh', 'w') as submit_bash:
    command = '''\
#! /bin/bash
# automatically generated submission file
# the defined standart parameters can be overwritten by command line parameters

MPI_PROCESS={}
VERBOSE={}
FRAMERATE={}
MESH_NAME={}
DTSTART={}
TEND={}
TIMEADAPTION=2
RICHARDSONSAFETY={}
RICHARDSONTOL={}
KRYLOVSOLVER={}
KRYLOVPREC={}
NEWTONABS={}
NEWTONREL={}

while [[ "$#" -gt 0 ]]; do case $1 in
  -n|--mpiprocess) MPI_PROCESS="$2"; shift;;
  -v|--verbose) VERBOSE="$2"; shift;;
  -fr|--framerate) FRAMERATE="$2"; shift;;
  # overwrite .config
  -m|--mesh) MESH_NAME="$2"; shift;;
  -dt|--dtstart) DTSTART="$2"; shift;;
  -T|--Tend) TEND="$2"; shift;;
  -type|--dttype) TIMEADAPTION="$2"; shift;;
  -rs|--richsafe) RICHARDSONSAFETY="$2"; shift;;
  -rt|--richtol) RICHARDSONTOL="$2"; shift;;
  -ls|--solver) KRYLOVSOLVER="$2"; shift;;
  -pc|--preconditioner) KRYLOVPREC="$2"; shift;;  
  -nr|--newtonrel) NEWTONREL="$2"; shift;;  
  -na|--newtonabs) NEWTONABS="$2"; shift;;  

  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# automatically added submissions
'''.format(mpiprocs, verbosity, framerate, mesh_name, dt, T_end, floatToString(rich_safe), floatToString(rich_tol),
           solver, preconditioner, floatToString(newton_abs),
           floatToString(newton_rel))
    submit_bash.write(command)
os.chmod(parent_path + '/submit-corner.sh', 0o755)

# fill the case directories with run files and populate sumbit-all script
fisher_path = "../../../bin/"
mesh_path = "../../../mesh"
runs_path = "output"
for c in case_dirs:
    # add to submission all file
    with open(parent_path + '/submit-all.sh', 'a') as submit_bash:
        command = '''{}/job.sh "$MESH_NAME" "$MPI_PROCESS" "$TEND" "$DTSTART" "$TIMEADAPTION" \\
        "$RICHARDSONSAFETY" "$RICHARDSONTOL" "$KRYLOVSOLVER" "$KRYLOVPREC" "$NEWTONABS" "$NEWTONREL" \\
        "$VERBOSE" "$FRAMERATE"
'''.format(c[0])
        submit_bash.write(command)

    # add to submission corner file if corner case
    if (c[1] == min(diffusion_pspace) or c[1] == max(diffusion_pspace)) and \
            (c[2] == min(rho_pspace) or c[2] == max(rho_pspace)) and \
            (c[3] == min(radius_pspace) or c[3] == max(radius_pspace)):
        print("corner case {}".format(c))
        with open(parent_path + '/submit-corner.sh', 'a') as submit_bash:
            command = '''{}/job.sh "$MESH_NAME" "$MPI_PROCESS" "$TEND" "$DTSTART" "$TIMEADAPTION" \\
                "$RICHARDSONSAFETY" "$RICHARDSONTOL" "$KRYLOVSOLVER" "$KRYLOVPREC" "$NEWTONABS" "$NEWTONREL" \\
                "$VERBOSE" "$FRAMERATE"
        '''.format(c[0])
            submit_bash.write(command)

    # write down job description
    with open(parent_path + "/" + c[0] + '/job.sh', 'w') as job_bash:
        description = '''
# this script is automatically generated by setupRun for {}
# job parameters:
#  path to executable = {}
#  date = {}
#  D_w = {}
#  rho = {}
#  radius = {}
#  diffusion_fac = {}  
#  diffusion_steps = {}  
#  rho_steps = {}  
#  radius_steps = {}
#  initial_condition = {}    
#
'''.format(c[0], fisher_path, date, c[1], c[2], c[3], diffusion_fac, diffusion_steps, rho_steps, radius_steps,
           initial_condition)

        command = '''\
#! /bin/bash
# run parameters
MESH_NAME=$1
MPI_PROCESS=$2
TEND=$3
DT_START=$4
RUNTYPE=$5
RICHSAFE=$6
RICHTOL=$7
SOLVER=$8
PREC=$9
NEWTONABS=${{10}}
NEWTONREL=${{11}}
VERBOSE=${{12}}
FRAMERATE=${{13}}

mpirun -n "$MPI_PROCESS" {}FisherSolver \\
    "{}"  "{}" "$MESH_NAME" "simulation-output" "out" "iterationdata" \\
    "{}" "{}" "{}" "{}" "1" \\
    "0.00000001" "$DT_START" "0.1" "$TEND" "$FRAMERATE" \\
    "{}" "{}" "{}" "1"\\
    "$VERBOSE" "$RUNTYPE" "$RICHTOL" "$RICHSAFE" "1"\\
    "$NEWTONABS" "$NEWTONREL" "-1" "-1" \\
    "$SOLVER" "$PREC" \\
    "0" "{}" "{}" "{}"
'''.format(fisher_path,
           c[0], mesh_path,
           initial_condition[0], initial_condition[1], initial_condition[2], c[3],
           c[1], c[1] / diffusion_fac, c[2],
           translation[0], translation[1], translation[2])
        # write file
        job_bash.write(description + command)

    # set execution permission
    os.chmod(parent_path + "/" + c[0] + '/job.sh', 0o755)
