# Dolfin Fisher Solver

### Installation requirements
Following packages have been used during developement of this software:

Ubuntu (18.04): OS  
GCC (6.3.0): Compiler  
CMake (minimum 3.5): Helps with compilation process   
FEniCS (2019.1): Automated solution of Differential Equations ([Installation guide](https://fenics.readthedocs.io/en/latest/installation.html))  
OpenMPI (3.0.1): Library for distributed computing, is integrated in FEniCS


### Build FisherSolver  
Follow these instructions to build FisherSolver executable:  
1. clone this repo to machine
2. make sure all requirements are installed and set-up
  * **Euler:** To build and run the simulation on the Euler cluster, load modules via the commands:  
    `env2lmod`  
    `module load gcc/6.3.0 cmake/3.15.3 openmpi/3.0.1 fenics/2019.1.0`
  * **local:** See the installation requirements
3. in DolfinFisherSolver parent folder, create missing directories:  
  `mkdir DolfinFisherSolver/build && mkdir output` 
4. Compile executable:
  * **Euler**: `cd DolfinFisherSolver && ./cmake_init.sh && cd build && make`
  * **local**: `cd DolfinFisherSolver/build && cmake .. && make`  
  
#### Euler cluster bug
Use the script `/DolfinFisherSolver/cmake_init.sh` to generate the CMake files, then use `make` as usual in the build folder.
The mesh partitioner of the Scotch library has a bug and is not working when the mesh reaches a certain number of elements, a fixed version of dolfin must be used. Source it with `source $HOME/usr/source_env.sh`   
If the correct version of dolfin is used can be check with:
`ldd FisherSolver | grep dolfin`

  
 ### Run FisherSolver
 After building, the executable should be in /bin.  
 The simulation is run via the script "runFisherSolver.sh". All the default simulation parameters are read-in from the "FisherSolver.config" file. 
If some of the parameters in the "FisherSolver.config" file are deleted or instantiated with unfeasible values, the behavior is undefined. There is the possibility to overwrite some values via command line parameters. The best is to take a quick look into "runFisherSolver.sh" to know what can be done.  
 
The provided /mesh folder is used as standart path for mesh files. This can be changed in "FisherSolver.config" and "runFisherSolver" files. 
Compatible formats are "h5" and "xml" , but it is recommended from the FEniCS community to use "h5".  
The "xml_to_h5.py" script can be used to transform .xml meshes to .h5 meshes. It takes the xml mesh name, for example "test-mesh.xml" as command line parameter and creates an .h5 version of the mesh with the name "test-mesh.h5"  
  

