# Dolfin Fisher Solver

### Installation and prerequisites
Following packages have been used during developement of this software:

Ubuntu (18.04): OS  
GCC (6.3.0): Compiler  
CMake (minimum 3.5): Helps with compilation process   
FEniCS (2019.1): Automated solution of Differential Equations ([Installation guide](https://fenics.readthedocs.io/en/latest/installation.html))
OpenMPI (3.0.1): Library for distributed computing, is integrated in FEniCS

#### ETH Euler cluster
To build and run the simulation on the Euler cluster, load following modules via the commands:  
`env2lmod`  
`module load gcc/6.3.0 cmake/3.15.3 openmpi/3.0.1 fenics/2019.1.0`  
Since the mesh partitioner of the Scotch library has a bug and is not working when the mesh reaches a certain number of elements, a fixed version of dolfin must be used. Source it with:  
`source /usr/source_env.sh`  
It then can be compiled with the correct version of dolfin. Use the script `/DolfinFisherSolver/cmake_init.sh` to generate the CMake files, then use `make` as usual in the build folder.
To make sure the correct version is used run commands:  
`ldd FisherSolver | grep dolfin`



### Build FisherSolver
1. clone this repo to local machine
2. make sure all requirements are installed
3. in DolfinFisherSolver parent folder run following commands:  
  `mkdir DolfinFisherSolver/build && mkdir output`  
  `cd DolfinFisherSolver/build && cmake .. && make`  
  
 ### Run FisherSolver
 After building, the executable should be in /bin.  
 The simulation is run via the script "runFisherSolver.sh". All the default simulation parameters are read-in from the "FisherSolver.config" file. The parameters in the "FisherSolver.config" file must be instantiated there. Otherwise some inputs to the FisherSolver are missing and the behavior is undefined. Try to work with feasible values, otherwise the behavior is again undefined. There is the possibility to overwrite some values via command line parameters. The best is to take a quick look into "runFisherSolver.sh" to know what can be done.  
 
 **Mesh**: The mesh defined in "FisherSolver.config" must be in the /mesh, folder. If not created externally, simple meshes can be created with the python script located in the /mesh folder.
 

 **Euler**: ` bsub -N -W 00:30 -R "rusage[scratch=10000, mem=12800]" -n 10 ./runFisherSolver.sh -n 10 -ofo "ten-cores"`  
 There can be used a maximum of 128000 MB of memory, therefore mem = 128000MB / n  

### Run weak scaling test
After building, a weak scaling test can be performed by running weakScalingTest.sh. This test submits a couple of batch jobs with processor numbers from 1 up to 48. The workload is once scaled in terms of timesteps, and once in terms of mesh resolution. Be sure to run "simple_mesh_creator.py" in advance to create the necessary meshes.

### Tensors
Each tensor is a dolfin::Expression. While solving a problem, the tensors eval()-function get's called and evaluated for every cell.
