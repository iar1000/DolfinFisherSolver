
# Performance Tester
This subpart of the simulation is testing different linear solver/ preconditioner pairs for the Newton solver. It is currently adapted to work on Euler.

### Build PerformanceTester
All the requirements from building FisherSolver must be fullfiled to build the PerformanceTester. If this is the case do the following:
1. build missing folders from within /Performance folder: 
 `mkdir src/build`
2. Setup and compile
`cd src && ./cmake_init.sh && cd build && make`

### Run PerformanceTester experimental test
For better understanding the problem that's being solved the type 4 test of the PerformanceTester can be run with the command line option `--type 4 --ls <desired solver> --pc <desired preconditioner>`.  
This allows to better inspect the underlying PETSc routines. Here is a variety of commands to log the krylov solver:
1. `--petsc.ksp_monitor`
2. `--petsc.ksp_monitor_true_residual`
3. `--petsc.ksp_monitor_singular_value`

### Run PerformanceTester weak scaling test
The runPerformanceTester.sh script runs the test for a range of processes.  
The script can be parameterised on the size of memory used and the number of dof's per core. If the number of dof's is to high, it can happen that the LSF memory limit is reached.  
The script produces xml and latex timetables with times of the different steps of the test.
