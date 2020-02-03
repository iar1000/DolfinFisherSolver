
#include <vector>
#include <dolfin.h>
#include <Eigen/Dense>

#include "TimeStepper.h"
#include "ReactionDiffusionProblem.h"
#include "Initializers.h"
#include "Tensors.h"
#include "ReaderWriter.h"


// create test value map for rect-100on100-res100.h test mesh
// @TODO ValueMapper
std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> get_101on101_test_cm()
{
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vm;

	// init valuemaps to form pattern
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vm_w = Eigen::MatrixXd::Ones(101, 101);
	vm_w.block<40, 101>(7, 0) = Eigen::MatrixXd::Zero(40, 101);
	vm_w.block<40, 101>(54, 0) = Eigen::MatrixXd::Zero(40, 101);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vm_g = Eigen::MatrixXd::Zero(101, 101);
	vm_g.block<40, 101>(7, 0) = Eigen::MatrixXd::Ones(40, 101);
	vm_g.block<40, 101>(54, 0) = Eigen::MatrixXd::Ones(40, 101);

	vm.push_back(vm_w);
	vm.push_back(vm_g);
	return vm;
};

int main(int argc, char* argv[]){

	// initialize MPI
	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// paths and names
	///////////////////////////////////
	std::string outputParent = argv[1];
	std::string meshParent = argv[2];
	std::string meshName = argv[3];
	std::string tagFolder = argv[4];
	std::string tagFile = argv[5];
	// initial condition
	////////////////////////////////////
	double cx = atof(argv[6]);
	double cy = atof(argv[7]);
	double cz = atof(argv[8]);
	double radius = atof(argv[9]);
	double value = atof(argv[10]);
	// timestepping parameters
	//////////////////////////////////
	double dt_min = atof(argv[11]);
	double dt = atof(argv[12]);
	double dt_max = atof(argv[13]);
	double T = atof(argv[14]);
	int framesPerTimeUnit = atoi(argv[15]);
	// reaction-diffusion coefficient
	////////////////////////////////////
	double Dw = atof(argv[16]);
	double Dg = atof(argv[17]);
	double rho = atof(argv[18]);
	// other
	///////////////////////////////////
	int verbose = atoi(argv[19]);
	int timeAdaption = atoi(argv[20]);


	// in/output handler
	ReaderWriter putput = ReaderWriter(rank, outputParent, meshParent);

	// mesh read-in
	std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>();
	std::pair<bool, int> meshInfo = putput.loadMesh(mesh, meshName);
	int dimensions = 0;
	if(meshInfo.first){	// check if mesh loaded successful
		dimensions = meshInfo.second;
	}
	else{
		return 0;
	}

	// create initial condition
	std::shared_ptr<dolfin::Expression> initialCondition;
	if(dimensions == 2){
		initialCondition = std::make_shared<InitializerCircle>(cx, cy, radius, value);
		// make dummy variable to print
		InitializerCircle dummy(cx, cy, radius, value);
		putput.addComponent(dummy.asString());
	}
	if(dimensions == 3){
		initialCondition = std::make_shared<InitializerSphere>(cx, cy, cz, radius, value);
		// make dummy variable to print
		InitializerSphere dummy(cx, cy, cz, radius, value);
		putput.addComponent(dummy.asString());
	}

	// create diffusion tensor
	std::shared_ptr<TensorConstant> DConstant = std::make_shared<TensorConstant>(rank, Dw);
	std::shared_ptr<TensorSpatial2D> DSpatial2D = std::make_shared<TensorSpatial2D>(rank, Dw, Dg, get_101on101_test_cm());
	putput.addComponent(DSpatial2D->asString());

	// create pde problem
	std::shared_ptr<ReactionDiffusionProblem> problem = std::make_shared<ReactionDiffusionProblem>(rank, mesh, DSpatial2D, rho, dt);
	putput.addComponent(problem->asString());

	// create solver
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>();
	solver->parameters["linear_solver"] = "lu";
	solver->parameters["convergence_criterion"] = "incremental";
	solver->parameters["maximum_iterations"] = 50;
	solver->parameters["relative_tolerance"] = 1e-10;
	solver->parameters["absolute_tolerance"] = 1e-10;

	// create time stepper
	TimeStepper timeStepper = TimeStepper(rank, problem, solver, T, dt_min, dt_max);
	putput.addComponent(timeStepper.asString());

	// create output file
	std::string fileType = "pvd";
	auto pathReturn = putput.getFilePath(tagFolder, tagFile, fileType);
	if(!pathReturn.first){	// check if path to file loaded successful
		return 0;
	}
	std::shared_ptr<dolfin::File> file = std::make_shared<dolfin::File>(pathReturn.second);

	// create pre-simulation info
	putput.createRunInfo(tagFolder, tagFile);

	// run simulation
	RuntimeTracker tracker3 = timeStepper.run(timeAdaption, verbose, initialCondition, file, framesPerTimeUnit, dt);

	// overwrite INFO file with post-simulation details
	putput.addComponent(tracker3.asString());
	putput.createRunInfo(tagFolder, tagFile);

	MPI_Finalize(); //seems to trigger an abort
}



