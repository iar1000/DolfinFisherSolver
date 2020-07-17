# Parameter Inference
## Setup
Generate a run directory and run scripty with `setupRun.py`.  
Take a look at the script to see the parameters you can set.  

The following are generated inside the run directory:
1. A folder for each case with a submission script `job.sh`.  
2. `submit-all.sh` to submit all cases to the lsf system.  
Parameters regarding the precisio, MPI processes and LSF execution time limit can be given as command line parameter, 
get detailed information by looking at the script.   
3. `submit-corner.sh` is analog to 2., but only submits corner cases.
4. `check-all.py` is printing an overview of the status of all cases  
5. `check-corner.py` is printing an overview of the status of the corner cases  
6. `clean-failed.py` is preparing the `submit-failed.sh` file and is cleaning the case directories of the stuck and timed-out
cases for a re-run. The lsf.out, lsf.err and all other files, including old simulation output is being **deleted**! So save the state you want to keep before running this script
Always run this script before re-submitting stuck cases.
7. `submit-failed.sh` is submiting the stuck and timed-out cases again. Consider changing the number of MPI processes to reduce the
probability of the case to get stuck again. (My experience)
