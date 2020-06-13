# script to generate a run enviroment for parametric inference
# automatically creates a run directory in /runs with generated cases

import sys
import os
from datetime import date
import numpy as np
import stat
import argparse


# https://stackoverflow.com/questions/2440692/formatting-floats-without-trailing-zeros
def floatToString(inputValue):
    return ('%.15f' % inputValue).rstrip('0').rstrip('.')

parser = argparse.ArgumentParser()
parser.add_argument("--Tend", type=float, help="T_end of the simulations", default=10)
parser.add_argument("--Dmin", type=float, help="minimum diffusion coefficient", default=0.1)
parser.add_argument("--Dmax", type=float, help="maximum diffusion coefficient", default=0.6)
parser.add_argument("--Dsteps", type=int, help="discretization steps between Dmin and Dmax", default=10)
parser.add_argument("--Rhomin", type=float, help="minimum diffusion coefficient", default=0.025)
parser.add_argument("--Rhomax", type=float, help="maximum diffusion coefficient", default=0.25)
parser.add_argument("--Rhosteps", type=int, help="discretization steps between Rhomin and Rhomax", default=10)
parser.add_argument("--Radiusmin", type=float, help="minimum initial radius", default=3)
parser.add_argument("--Radiusmax", type=float, help="maximum initial radius", default=3)
parser.add_argument("--Radiussteps", type=int, help="discretization steps between Radiusmin and Radiusmax", default=1)

parser.add_argument("--procs", type=int, help="number mpi processors (must be multiple of 24/36)", default=24)
parser.add_argument("--runtime", type=str, help="lsf time limit", default="04:00")
parser.add_argument("--mesh", type=str, help="name of the mesh", default="lh-plial-dof-5500k.h5")

args = parser.parse_args()

# Parameter inference
diffusion_fac = 10  # factor D_w / D_g
diffusion_min = args.Dmin           # min D_W
diffusion_max = args.Dmax           # max D_W
diffusion_steps = args.Dsteps       # discretization steps of parameter range
rho_min = args.Rhomin               # min rho
rho_max = args.Rhomax               # max rho
rho_steps = args.Rhosteps           # discretization steps of parameter range
radius_min = args.Radiusmin         # min r0
radius_max = args.Radiusmax         # max r0
radius_steps = args.Radiussteps     # discretization steps of parameter range

# Fixed Arguments
# Runtime arguments
mpiprocs = args.procs
lsftime = args.runtime
verbosity = 3
# Model arguments
date = date.today().strftime("%d-%m-%Y")
initial_condition = [46, 132, 67, 3]  # x, y, z, radius
translation = [17, 21, 30]
# Solver arguments
mesh_name = args.mesh
T_end = args.Tend
dt = 0.0001
solver = "cg"
preconditioner = "hypre_euclid"
rich_tol = 0.000001
rich_safe = 0.95
newton_abs = 0.00000001
newton_rel = 0.00000001
framerate = 1

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
radius_stepsize = (1 if radius_range == 0 else radius_range / radius_steps)
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
corner_dirs = []
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
TIME={}

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
  -t|--lsftime) TIME="$2"; shift;; 
  
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# automatically added submissions
'''.format(mpiprocs, verbosity, framerate, mesh_name, dt, T_end, floatToString(rich_safe), floatToString(rich_tol),
           solver, preconditioner, floatToString(newton_abs),
           floatToString(newton_rel), lsftime)
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
TIME={}

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
  -t|--lsftime) TIME="$2"; shift;;   

  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# automatically added submissions
'''.format(mpiprocs, verbosity, framerate, mesh_name, dt, T_end, floatToString(rich_safe), floatToString(rich_tol),
           solver, preconditioner, floatToString(newton_abs),
           floatToString(newton_rel), lsftime)
    submit_bash.write(command)
os.chmod(parent_path + '/submit-corner.sh', 0o755)

# create submission file for corner cases
with open(parent_path + '/submit-failed.sh', 'w') as submit_bash:
    command = '''\
#! /bin/bash
# automatically generated submission file for failed cases
# change the number of processors
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
TIME={}

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
  -t|--lsftime) TIME="$2"; shift;;   

  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# automatically added submissions by check-all.py
# SENTINEL LINE
'''.format(mpiprocs, verbosity, framerate, mesh_name, dt, T_end, floatToString(rich_safe), floatToString(rich_tol),
           solver, preconditioner, floatToString(newton_abs),
           floatToString(newton_rel), lsftime)
    submit_bash.write(command)
os.chmod(parent_path + '/submit-failed.sh', 0o755)


# fill the case directories with run files and populate sumbit-all script
fisher_path = "../../../bin/"
mesh_path = "../../../mesh"
lsf_output_name = "lsf-output"
for c in case_dirs:
    # add to submission all file
    with open(parent_path + '/submit-all.sh', 'a') as submit_bash:
        command = '''{}/job.sh "$MESH_NAME" "$MPI_PROCESS" "$TEND" "$DTSTART" "$TIMEADAPTION" \\
        "$RICHARDSONSAFETY" "$RICHARDSONTOL" "$KRYLOVSOLVER" "$KRYLOVPREC" "$NEWTONABS" "$NEWTONREL" \\
        "$VERBOSE" "$FRAMERATE" "$TIME"
'''.format(c[0])
        submit_bash.write(command)

    # add to submission corner file if corner case
    if (c[1] == min(diffusion_pspace) or c[1] == max(diffusion_pspace)) and \
            (c[2] == min(rho_pspace) or c[2] == max(rho_pspace)) and \
            (c[3] == min(radius_pspace) or c[3] == max(radius_pspace)):
        print("corner case {}".format(c))
        corner_dirs.append(c)
        with open(parent_path + '/submit-corner.sh', 'a') as submit_bash:
            command = '''{}/job.sh "$MESH_NAME" "$MPI_PROCESS" "$TEND" "$DTSTART" "$TIMEADAPTION" \\
                "$RICHARDSONSAFETY" "$RICHARDSONTOL" "$KRYLOVSOLVER" "$KRYLOVPREC" "$NEWTONABS" "$NEWTONREL" \\
                "$VERBOSE" "$FRAMERATE" "$TIME"
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
# LSF parameters
MPI_PROCESS=$2
TIME=${{14}}


 bsub -N -W "$TIME" -n "$MPI_PROCESS" -o {}/lsf.out -e {}/lsf.err -R fullnode mpirun -n "$MPI_PROCESS" {}FisherSolver \\
    "{}"  "{}" "$MESH_NAME" "simulation-output" "out" "iterationdata" \\
    "{}" "{}" "{}" "{}" "1" \\
    "0.00000001" "$DT_START" "0.1" "$TEND" "$FRAMERATE" \\
    "{}" "{}" "{}" "1"\\
    "$VERBOSE" "$RUNTYPE" "$RICHTOL" "$RICHSAFE" "1"\\
    "$NEWTONABS" "$NEWTONREL" "-1" "-1" \\
    "$SOLVER" "$PREC" \\
    "0" "{}" "{}" "{}"
'''.format(c[0], c[0], fisher_path,
           c[0], mesh_path,
           initial_condition[0], initial_condition[1], initial_condition[2], c[3],
           c[1], c[1] / diffusion_fac, c[2],
           translation[0], translation[1], translation[2])
        # write file
        job_bash.write(description + command)

    # set execution permission
    os.chmod(parent_path + "/" + c[0] + '/job.sh', 0o755)
    #

###############################################################################
####################### CHECK AND RERUN #######################################
###############################################################################

# create check status of corner cases
corner_dirs_names = [x[0] for x in corner_dirs]
with open(parent_path + '/check-corner.py', 'w') as submit_bash:
    command = '''\
# automatically generated script to check status of corner case submissions

import os
from datetime import datetime

corner_directories = {}
all_directories = os.listdir()

# collect status infos
pending_cases = []  # cases no submission has been made, yet
done_cases = []  # cases at least one sucessfull simulation has been run
midrun_cases = []  # cases that are mid-simulation
init_cases = []  # case that are currently initalizing the simulation (mesh read in)
timedout_cases = []  # cases that are lsf timed out
stuck_cases = []  # cases that got stuck in mesh read-in
failed_cases = [] # cases that got submitted but failed without starting the simulation
# check all cases
for dir in all_directories:
    if dir in corner_directories:

        has_done = False  # indicator if at least one run has been completed
        has_started = True  # indicator if at least one run has started
        case_statuses = []  # collect data of all runs of this case

        path = dir + "/simulation-output/"
        # check if simulation output folder exists, created by FisherSolver
        # get all files in this folder
        lsf_directory = []
        try:
            case_directories = os.listdir(path)
            lsf_directory = os.listdir(dir)
        except:
            case_directories = []
            # check if it wasn't submitted or simply failed
            lsfoutput = [x for x in lsf_directory if ("lsf.out" in x)]
            if lsfoutput:
                failed_cases.append([dir])
            else:
                pending_cases.append([dir])
            has_started = False

        # catch pending cases with no submissions yet
        if not has_started:
            continue

        # collect status of different processor runs
        case_status_files = [x for x in case_directories if ("_STATUS-" in x)]
        for sf in case_status_files:
            nprocs = sf.split("-")[1].split(".")[0]
            # check status files
            with open(path + sf) as file:
                lines = file.readlines()
                if lines:
                    start_time = datetime.strptime(lines[1].replace("\\x00", "").replace("\\n", "")[-19:],
                                                   '%Y-%m-%d %H:%M:%S')
                    status = " ".join(lines[-1].split())
                    status_name = status.split("-")[0]
                    status_time = datetime.strptime(lines[-1].replace("\\x00", "").replace("\\n", "")[-19:],
                                                    '%Y-%m-%d %H:%M:%S')
                    case_statuses.append([nprocs, start_time, status_name, status_time])
                else:
                    has_started = False #@TODO

        # check if at least one simulation was sucessfull
        for run_status in case_statuses:
            if "Simulation 4 done" in run_status[2]:
                has_done = True
                done_cases.append([dir, run_status[0], abs((run_status[3] - run_status[1])), len(case_statuses)])
                break
        if has_done:
            continue

        # at this point all pending or sucessefull finnished cases should be handled and in the corresponding folders
        # check if the simulation for this directory is currently running or failed
        lsfoutput = [x for x in lsf_directory if ("lsf.out" in x)]
        case_statuses = sorted(case_statuses, reverse=True, key=lambda x: x[3])
        
        # failed becasue not case status is available but there is a simulation-output folder
        if not case_statuses:
            failed_cases.append([dir])
            continue
            
        # lsf output in directory means case has failed
        if lsfoutput:
            if "Simulation 3 run" in case_statuses[0][2]:
                timedout_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            elif "Simulation 2 mesh" in case_statuses[0][2]:
                stuck_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            else:
                print("what happend here? ", case_statuses[0][2])
        else:
            if "Simulation 3 run" in case_statuses[0][2]:
                midrun_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            elif "Simulation 2 mesh" in case_statuses[0][2]:
                init_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            else:
                print("what happend here? ", case_statuses[0][2])



# print status infos
print("Status infos:")
print("\\tDone cases : {{}}\\n"
      "\\tFailed cases: {{}}\\n"
      "\\tInit cases : {{}}\\n"
      "\\tRunning cases : {{}}\\n"
      "\\tTimed-out cases : {{}}\\n"
      "\\tStuck cases : {{}}\\n"
      "\\tPending cases : {{}}".format(len(done_cases), len(failed_cases), len(init_cases), len(midrun_cases),
                                    len(timedout_cases), len(stuck_cases), len(pending_cases)))

# failed cases
print("\\nFailed cases:")
print('%-35s' % ("case directory"))
print("--------------------------------------------------------")
failed_cases = sorted(failed_cases, key=lambda x: x)
for i in range(len(failed_cases)):
    print('%-35s' % (failed_cases[i][0]))

# stuck cases
print("\\nStuck cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
stuck_cases = sorted(stuck_cases, key=lambda x: x)
for i in range(len(stuck_cases)):
    print('%-35s %-5i %-5i' % (stuck_cases[i][0], int(stuck_cases[i][1]), int(stuck_cases[i][2])))

# timed out cases
print("\\nTimed-out cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
timedout_cases = sorted(timedout_cases, key=lambda x: (x[0], x[1]))
for i in range(len(timedout_cases)):
    print('%-35s %-5i %-5i' % (timedout_cases[i][0], int(timedout_cases[i][1]), int(timedout_cases[i][2])))

# init cases
print("\\nInit cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
init_cases = sorted(init_cases, key=lambda x: (x[0], x[1]))
for i in range(len(init_cases)):
    print('%-35s %-5i %-5i' % (init_cases[i][0], int(init_cases[i][1]), int(init_cases[i][2])))

# running cases
print("\\nRunning cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
midrun_cases = sorted(midrun_cases, key=lambda x: (x[0], x[1]))
for i in range(len(midrun_cases)):
    print('%-35s %-5i %-5i' % (midrun_cases[i][0], int(midrun_cases[i][1]), int(midrun_cases[i][2])))

# pending cases
print("\\nPending cases:")
print('%-35s' % ("case directory"))
print("--------------------------------------------------------")
pending_cases = sorted(pending_cases, key=lambda x: x)
for i in range(len(pending_cases)):
    print('%-35s' % (pending_cases[i][0]))

# done cases
print("\\nDone cases:")
print('%-35s %-6s %-10s %-6s' % ("case directory", "procs", "runtime", "tries"))
print("--------------------------------------------------------")
done_cases = sorted(done_cases, key=lambda x: (x[0], x[1]))
for i in range(len(done_cases)):
    print('%-35s %-5i %-10s %-5i' % (done_cases[i][0], int(done_cases[i][1]), done_cases[i][2], int(done_cases[i][3])))
'''.format(corner_dirs_names)
    submit_bash.write(command)
os.chmod(parent_path + '/check-corner.py', 0o755)

# create check status of all cases
case_dirs_names = [x[0] for x in case_dirs]
with open(parent_path + '/check-all.py', 'w') as submit_bash:
    command = '''\
# automatically generated script to check status of all case submissions

import os
from datetime import datetime

corner_directories = {}
all_directories = os.listdir()

# collect status infos
pending_cases = []  # cases no submission has been made, yet
done_cases = []  # cases at least one sucessfull simulation has been run
midrun_cases = []  # cases that are mid-simulation
init_cases = []  # case that are currently initalizing the simulation (mesh read in)
timedout_cases = []  # cases that are lsf timed out
stuck_cases = []  # cases that got stuck in mesh read-in
failed_cases = [] # cases that got submitted but failed without starting the simulation
# check all cases
for dir in all_directories:
    if dir in corner_directories:

        has_done = False  # indicator if at least one run has been completed
        has_started = True  # indicator if at least one run has started
        case_statuses = []  # collect data of all runs of this case

        path = dir + "/simulation-output/"
        # check if simulation output folder exists, created by FisherSolver
        # get all files in this folder
        lsf_directory = []
        try:
            case_directories = os.listdir(path)
            lsf_directory = os.listdir(dir)
        except:
           # check if it wasn't submitted or simply failed
            lsfoutput = [x for x in lsf_directory if ("lsf.out" in x)]
            if lsfoutput:
                failed_cases.append([dir])
            else:
                pending_cases.append([dir])
            has_started = False

        # catch pending cases with no submissions yet
        if not has_started:
            continue

        # collect status of different processor runs
        case_status_files = [x for x in case_directories if ("_STATUS-" in x)]
        for sf in case_status_files:
            nprocs = sf.split("-")[1].split(".")[0]
            # check status files
            with open(path + sf) as file:
                lines = file.readlines()
                if lines:
                    start_time = datetime.strptime(lines[1].replace("\\x00", "").replace("\\n", "")[-19:],
                                                   '%Y-%m-%d %H:%M:%S')
                    status = " ".join(lines[-1].split())
                    status_name = status.split("-")[0]
                    status_time = datetime.strptime(lines[-1].replace("\\x00", "").replace("\\n", "")[-19:],
                                                    '%Y-%m-%d %H:%M:%S')
                    case_statuses.append([nprocs, start_time, status_name, status_time])
                else:
                    has_started = False #@TODO

        # check if at least one simulation was sucessfull
        for run_status in case_statuses:
            if "Simulation 4 done" in run_status[2]:
                has_done = True
                done_cases.append([dir, run_status[0], abs((run_status[3] - run_status[1])), len(case_statuses)])
                break
        if has_done:
            continue

        # at this point all pending or sucessefull finnished cases should be handled and in the corresponding folders
        # check if the simulation for this directory is currently running or failed
        lsfoutput = [x for x in lsf_directory if ("lsf.out" in x)]
        case_statuses = sorted(case_statuses, reverse=True, key=lambda x: x[3])
        
        # failed becasue not case status is available but there is a simulation-output folder
        if not case_statuses:
            failed_cases.append([dir])
            continue
            
        # lsf output in directory means case has failed
        if lsfoutput:
            if "Simulation 3 run" in case_statuses[0][2]:
                timedout_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            elif "Simulation 2 mesh" in case_statuses[0][2]:
                stuck_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            else:
                print("what happend here? ", case_statuses[0][2])
        else:
            if "Simulation 3 run" in case_statuses[0][2]:
                midrun_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            elif "Simulation 2 mesh" in case_statuses[0][2]:
                init_cases.append([dir, case_statuses[0][0], len(case_statuses)])
            else:
                print("what happend here? ", case_statuses[0][2])



# print status infos
print("Status infos:")
print("\\tDone cases : {{}}\\n"
      "\\tFailed cases: {{}}\\n"
      "\\tInit cases : {{}}\\n"
      "\\tRunning cases : {{}}\\n"
      "\\tTimed-out cases : {{}}\\n"
      "\\tStuck cases : {{}}\\n"
      "\\tPending cases : {{}}".format(len(done_cases), len(failed_cases), len(init_cases), len(midrun_cases),
                                    len(timedout_cases), len(stuck_cases), len(pending_cases)))
# failed cases
print("\\nFailed cases:")
print('%-35s' % ("case directory"))
print("--------------------------------------------------------")
failed_cases = sorted(failed_cases, key=lambda x: x)
for i in range(len(failed_cases)):
    print('%-35s' % (failed_cases[i][0]))

# stuck cases
print("\\nStuck cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
stuck_cases = sorted(stuck_cases, key=lambda x: x)
for i in range(len(stuck_cases)):
    print('%-35s %-5i %-5i' % (stuck_cases[i][0], int(stuck_cases[i][1]), int(stuck_cases[i][2])))

# timed out cases
print("\\nTimed-out cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
timedout_cases = sorted(timedout_cases, key=lambda x: (x[0], x[1]))
for i in range(len(timedout_cases)):
    print('%-35s %-5i %-5i' % (timedout_cases[i][0], int(timedout_cases[i][1]), int(timedout_cases[i][2])))

# init cases
print("\\nInit cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
init_cases = sorted(init_cases, key=lambda x: (x[0], x[1]))
for i in range(len(init_cases)):
    print('%-35s %-5i %-5i' % (init_cases[i][0], int(init_cases[i][1]), int(init_cases[i][2])))

# running cases
print("\\nRunning cases:")
print('%-35s %-6s %-6s' % ("case directory", "procs", "tries"))
print("--------------------------------------------------------")
midrun_cases = sorted(midrun_cases, key=lambda x: (x[0], x[1]))
for i in range(len(midrun_cases)):
    print('%-35s %-5i %-5i' % (midrun_cases[i][0], int(midrun_cases[i][1]), int(midrun_cases[i][2])))

# pending cases
print("\\nPending cases:")
print('%-35s' % ("case directory"))
print("--------------------------------------------------------")
pending_cases = sorted(pending_cases, key=lambda x: x)
for i in range(len(pending_cases)):
    print('%-35s' % (pending_cases[i][0]))

# done cases
print("\\nDone cases:")
print('%-35s %-6s %-10s %-6s' % ("case directory", "procs", "runtime", "tries"))
print("--------------------------------------------------------")
done_cases = sorted(done_cases, key=lambda x: (x[0], x[1]))
for i in range(len(done_cases)):
    print('%-35s %-5i %-10s %-5i' % (done_cases[i][0], int(done_cases[i][1]), done_cases[i][2], int(done_cases[i][3])))
'''.format(case_dirs_names)
    submit_bash.write(command)
os.chmod(parent_path + '/check-all.py', 0o755)

# script to set up submit-stuck.sh submission script for stuck cases
# create check status of all cases
case_dirs_names = [x[0] for x in case_dirs]
with open(parent_path + '/clean-failed.py', 'w') as submit_bash:
    command = '''\
# automatically generated script to setup the submit-failed.sh script with job submissions of failed cases
# deletes the lsf output
# deletes simulation-output files (_STATUS-* files are kept for information about past runs)

import os
from datetime import datetime
import shutil

corner_directories = {}

# collect status infos
pending_cases = []  # cases that are waiting to be submitted or runnning
done_cases = []  # cases at least one sucessfull simulation has been run
failed_cases = []  # cases that have not sucessfully finished the simulation in the given time or failed

# check all cases
for dir in corner_directories:
    if dir in corner_directories:

        has_done = False  # indicator if at least one run has been completed
        has_started = True  # indicator if at least one run has started
        case_statuses = []  # collect data of all runs of this case

        # check if simulation output folder exists, created by FisherSolver
        path = dir + "/simulation-output/"
        lsf_directory = []
        try:
            case_directories = os.listdir(path)
        except:
            lsf_directory = os.listdir(dir)
            # check if it wasn't submitted or simply failed
            lsfoutput = [x for x in lsf_directory if ("lsf.out" in x)]
            if lsfoutput:
                failed_cases.append([dir])          # FAILED
            else:
                pending_cases.append([dir])         # PENDING 1. SUBMIT
            continue

        # HAS SIMULATION-OUTPUT FOLDER

        # collect status of different processor runs
        case_status_files = [x for x in case_directories if ("_STATUS-" in x)]
        for sf in case_status_files:
            nprocs = sf.split("-")[1].split(".")[0]
            # check status files
            with open(path + sf) as file:
                lines = file.readlines()
                if lines:
                    start_time = datetime.strptime(lines[1].replace("\\x00", "").replace("\\n", "")[-19:],
                                                   '%Y-%m-%d %H:%M:%S')
                    status = " ".join(lines[-1].split())
                    status_name = status.split("-")[0]
                    status_time = datetime.strptime(lines[-1].replace("\\x00", "").replace("\\n", "")[-19:],
                                                    '%Y-%m-%d %H:%M:%S')
                    case_statuses.append([nprocs, start_time, status_name, status_time])
                else:
                    has_started = False

        # check if at least one simulation was sucessfull
        for run_status in case_statuses:
            if "Simulation 4 done" in run_status[2]:
                has_done = True
                done_cases.append([dir, run_status[0], abs((run_status[3] - run_status[1])), len(case_statuses)])
                break
        if has_done:
            continue

        # ALL CASES THAT DID NOT FINISHED AT LEAST ONE SIMULATION SUCESSFULLY

        # check if the simulation for this directory is currently running or failed
        lsfoutput = [x for x in lsf_directory if ("lsf.out" in x)]
        # case was either stuck, run out of time or failed
        if lsfoutput:
            failed_cases.append([dir])
        # case is pending
        else:
            pending_cases.append([dir])  # PENDING NOT 1. SUBMIT (Has simulation-out folder but didn't finish)

# print status infos
print("Setup re-run of failed cases:")
print("\\tFailed cases : {{}}  ".format(len(failed_cases)))
print("\\tPending cases : {{}}  ".format(len(pending_cases)))
print("\\tDone cases : {{}}  ".format(len(done_cases)))


# remove old simulation data
for case in failed_cases:
    try:
        path = case[0] + "/simulation-output/"
        os.listdir(path)
    except:
        print("WARNING: should not be here, all simulations w/o simulation-output folder should be in Pending")
        continue
    # Remove old pvd files
    try:
        shutil.rmtree(path + "pvd")
    except:
        print("ERROR: failed removing /pvd of ", case[0])
    # remove old function files
    try:
        shutil.rmtree(path + "functions")
    except:
        print("ERROR: failed removing /functions of ", case[0])
    # remove INFO file
    try:
        os.remove(path + "_INFO-out.txt")
    except:
        print("ERROR: failed removing _INFO of ", case[0])
    # remove iteration
    try:
        os.remove(path + "iterationdata.csv-0")
    except:
        print("ERROR: failed removing iteration.csv of ", case[0])
    # remove timings
    try:
        os.remove(path + "timings-latex.txt")
    except:
        print("ERROR: failed removing timings of ", case[0])
    # remove lsf.out
    try:
        os.remove(path + "lsf.out")
    except:
        print("ERROR: failed removing lsf.out of ", case[0])
    # remove lsf.err
    try:
        os.remove(path + "lsf.err")
    except:
        print("ERROR: failed removing lsf.err of ", case[0])

# add cases to submit script
with open('submit-failed.sh', 'r+') as submit_bash:
    # read in old content of the file and clear afterwards
    old_file = submit_bash.read()
    submit_bash.seek(0)
    submit_bash.truncate()
    # delete all old submission jobs
    sep = "# SENTINEL LINE"
    preface = old_file.split(sep)[0] + sep + "\\n"
    submit_bash.write(preface)
    # append stuck submission jobs
    for case in failed_cases:
        command = \'\'\'
{{}}/job.sh "$MESH_NAME" "$MPI_PROCESS" "$TEND" "$DTSTART" "$TIMEADAPTION" \\
       "$RICHARDSONSAFETY" "$RICHARDSONTOL" "$KRYLOVSOLVER" "$KRYLOVPREC" "$NEWTONABS" "$NEWTONREL" \\
        "$VERBOSE" "$FRAMERATE" "$TIME"
        \'\'\'.format(case[0])
        submit_bash.write(command)
'''.format(case_dirs_names)
    submit_bash.write(command)
os.chmod(parent_path + '/clean-failed.py', 0o755)
