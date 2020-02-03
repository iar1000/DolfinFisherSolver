# Dolfin Fisher Solver

### Installation and prerequisites
Following packages have been used during developement of this software:

Ubuntu (18.04): OS  
CMake (minimum 3.5): Compilation process   
FEniCS (2019.1): Automated solution of Differential Equations ([Installation guide](https://fenics.readthedocs.io/en/latest/installation.html))

#### ETH Euler cluster
To build the simulation on the Euler cluster, load following modules via the commands:  
`env2lmod`  
`module load gcc/6.3.0 cmake/3.15.3 openmpi/3.0.1 fenics/2019.1.0`


### Build FisherSolver
1. clone this repo to local machine
2. in DolfinFisherSolver parent folder run following commands:  
  `mkdir DolfinFisherSolver/build && mkdir output`  
  `cd DolfinFisherSolver/build && cmake .. && make`  
  
 ### Run FisherSolver
 After building, the executable should be in /bin.  
 The simulation is run via the script "runFisherSolver.sh". All the default simulation parameters are read-in from the "FisherSolver.config" file. Ther is the possibility to overwrite some values via command line parameters. The best is to take a quick look into "runFisherSolver.sh" to know what can be done.  
 
 **Mesh**: The mesh defined in "FisherSolver.config" must be in the /mesh, folder. If not created externally, simple meshes can be created with the python script located in the /mesh folder.
 
**Note**: The parameters in the "FisherSolver.config" file must be instantiated there. Otherwise some inputs are missing and the behavior is undefined. Try to work with feasible values, otherwise the behavior is again undefined

 **Euler**: `bsub -N -W 24:00 -n 48 -R fullnode < runFisherSolver.sh -n 48 -ofo "not-default-name"`


### Tensors
Each tensor is a dolfin::Expression. While solving a problem, the tensors eval()-function get's called and evaluated for every cell.
