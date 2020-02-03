
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

	ReaderWriter putput = ReaderWriter(rank, "../output", "../mesh");

	// mesh
	std::string meshName = "rect-100on100-res-100.h5";
	std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>();
	std::pair<bool, int> meshInfo = putput.loadMesh(mesh, meshName);
	int dimensions = 0;
	if(meshInfo.first){
		dimensions = meshInfo.second;
	}
	else{
		return 0;
	}

	// timestepping parameters
	double dt_min = 0.0000001;
	double dt = 0.001;
	double dt_max = 0.001;
	double T = 0.1;

	// initial condition
	std::shared_ptr<dolfin::Expression> initialCondition;
	double cx = 50;
	double cy = 50;
	double cz = 5;
	double radius = 4;
	double value = 1;
	if(dimensions == 2){
		initialCondition = std::make_shared<InitializerCircle>(cx, cy, radius, value);
	}
	if(dimensions == 3){
		initialCondition = std::make_shared<InitializerSphere>(cx, cy, cz, radius, value);
	}
	//putput.addComponent(initialCondition->asString());

	// reaction coefficient
	double rho = 0.025;

	// diffusion tensor
	double Dw = 0.013;
	double Dg = 0.0013;
	std::shared_ptr<TensorConstant> DConstant = std::make_shared<TensorConstant>(rank, Dw);
	std::shared_ptr<TensorSpatial2D> DSpatial2D = std::make_shared<TensorSpatial2D>(rank, Dw, Dg, get_101on101_test_cm());
	putput.addComponent(DSpatial2D->asString());

	// pde problem
	std::shared_ptr<ReactionDiffusionProblem> problem = std::make_shared<ReactionDiffusionProblem>(rank, mesh, DSpatial2D, rho, dt);
	putput.addComponent(problem->asString());

	// solver
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>();
	solver->parameters["linear_solver"] = "lu";
	solver->parameters["convergence_criterion"] = "incremental";
	solver->parameters["maximum_iterations"] = 50;
	solver->parameters["relative_tolerance"] = 1e-10;
	solver->parameters["absolute_tolerance"] = 1e-10;

	// output parameters
	std::string tag_folder = "test-test";
	std::string tag_file = "out";
	std::string file_datatype = "pvd";
	auto pathReturn = putput.getFilePath(tag_folder, tag_file, file_datatype);
	if(!pathReturn.first){
		return 0;
	}
	std::shared_ptr<dolfin::File> file = std::make_shared<dolfin::File>(pathReturn.second);

	int framesPerTimeUnit = 60;

	// time stepper
	TimeStepper timeStepper = TimeStepper(rank, problem, solver, T, dt_min, dt_max);
	putput.addComponent(timeStepper.asString());

	// create pre-simulation info and run simulation
	putput.createRunInfo(tag_folder, tag_file);
	RuntimeTracker tracker3 = timeStepper.run(1, 3, initialCondition, file, framesPerTimeUnit, dt);

	// overwrite INFO file with post-simulation details
	putput.addComponent(tracker3.asString());
	putput.createRunInfo(tag_folder, tag_file);

	MPI_Finalize(); //seems to trigger an abort
}



