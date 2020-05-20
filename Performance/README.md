
# Performance Tester
This subpart of the simulation is testing different linear solver/ preconditioner pairs for the Newton solver. It is currently adapted to work on Euler.

### Build PerformanceTester
All the requirements from building FisherSolver must be fullfiled to build the PerformanceTester. If this is the case do the following:
1. build missing folders from within /Performance folder: 
 `mkdir src/build && mkdir output`
2. Setup and compile
`cd src && ./cmake_init.sh && cd build && make`

### Run PerformanceTester experimental test
For better understanding the problem that's being solved the type 4 test of the PerformanceTester can be run with the command line option `--type 1 --ls <desired solver> --pc <desired preconditioner> --meshname <name of mesh in ../../mesh>`.  
To better inspect the underlying PETSc routines check the PETSc [manual](https://www.mcs.anl.gov/petsc/documentation/index.html). Here are some commands to log the krylov solver:
1. `--petsc.ksp_monitor`
2. `--petsc.ksp_monitor_true_residual`
3. `--petsc.ksp_monitor_singular_value`

### Run scripts
The run scripts submit lsf jobs of PerformanceTester with different parameters.  
The names are self explaining. I used them on the Euler cluster, no guarantees given that they work anywhere else. The data generated lands in the /output folder.  
The generated data can be plotted with the scripts in /plotting. These scripts look for performance data files with the appropriate name in the /output folder. Only plot-iterationdata.py is looking for a file called "iterationdata.csv-0" in the local folder.
