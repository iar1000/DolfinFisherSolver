
# Performance Tester
This subpart of the simulation is testing different linear solver/ preconditioner pairs for the Newton solver. It is currently adapted to work on Euler.

### Build PerformanceTester
All the requirements from building FisherSolver must be fullfiled to build the PerformanceTester. If this is the case do the following in the Performance/ folder:
1. Setup and compile normal mode
`cd src && ./cmake_init.sh && cd build && make` or  
Setup and compile debug mode with gprof output
`cd src && ./cmake_init_debug.sh && cd build && make`  

### Run PerformanceTester experimental test
To get a performance overview of the available krylov space solver and preconditioner pairs, the type 1 test of the PerformanceTester can be run with the command line option `--type 1 --meshname <name of mesh in ../../mesh>`.  
The type 2 test is specified to the production enviroment and is only testing the Conjugate Gradient method: `--type 2 --meshname <name of mesh in ../../mesh>`   

To better inspect the underlying PETSc routines check the PETSc [manual](https://www.mcs.anl.gov/petsc/documentation/index.html). 
To get output from the PETSC library, the code must be setup and compiled with the non-debug `cmake_init.sh` script.  
Here are some commands to log the krylov solver:
1. `--petsc.ksp_monitor`
2. `--petsc.ksp_monitor_true_residual`
3. `--petsc.ksp_monitor_singular_value`

To get the PETSc profiling results, use `--petsc.log_view > petsc-profile.txt`

### Run scripts
The run scripts submit lsf jobs of PerformanceTester with different parameters.  
The names are self explaining. I used them on the Euler cluster, no guarantees given that they work anywhere else. The data generated lands in the /output folder.  
The generated data can be plotted with the scripts in /plotting. These scripts look for performance data files with the appropriate name in the /output folder. Only plot-iterationdata.py is looking for a file called "iterationdata.csv-0" in the local folder.
