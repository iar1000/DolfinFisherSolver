
#include <vector>
#include <dolfin.h>
#include <Eigen/Dense>

#include "TimeStepper.h"
#include "ReactionDiffusionProblem.h"
#include "Initializers.h"
#include "Tensors.h"
#include "ReaderWriter.h"

// create test value map for rect-10on10
std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> get_10on10_test_cm()
{
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vm;

	// init valuemaps to form pattern
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vm_w = Eigen::MatrixXd::Ones(11, 11);
	vm_w.block<5, 11>(0, 0) = Eigen::MatrixXd::Zero(5, 11);

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vm_g = Eigen::MatrixXd::Zero(11, 11);
	vm_g.block<5, 11>(0, 0) = Eigen::MatrixXd::Ones(5, 11);

	vm.push_back(vm_w);
	vm.push_back(vm_g);
	return vm;
};

// create test value map for box-10on10on10-res15.h test mesh
std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> get_10on10on10_test_cm(){
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vm_w;
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vm_g;

	for(int i = 0; i < 16; i++){
		auto ok = get_10on10_test_cm();
		vm_w.push_back(ok[0]);
		vm_g.push_back(ok[1]);
	}

	std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> vm;
	vm.push_back(vm_w);
	vm.push_back(vm_g);

	return vm;
}


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
	std::string tagCsv = argv[6];
	// initial condition
	////////////////////////////////////
	double cx = atof(argv[7]);
	double cy = atof(argv[8]);
	double cz = atof(argv[9]);
	double radius = atof(argv[10]);
	double value = atof(argv[11]);
	// timestepping parameters
	//////////////////////////////////
	double dt_min = atof(argv[12]);
	double dt = atof(argv[13]);
	double dt_max = atof(argv[14]);
	double T = atof(argv[15]);
	int framesPerTimeUnit = atoi(argv[16]);
	// reaction-diffusion coefficient
	////////////////////////////////////
	double Dw = atof(argv[17]);
	double Dg = atof(argv[18]);
	double rho = atof(argv[19]);
	double theta = atof(argv[20]);
	// other
	///////////////////////////////////
	int verbose = atoi(argv[21]);
	int timeAdaption = atoi(argv[22]);
	double richTol = atof(argv[23]);
	double richSafe = atof(argv[24]);



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
	else if(dimensions == 3){
		initialCondition = std::make_shared<InitializerSphere>(cx, cy, cz, radius, value);
		// make dummy variable to print
		InitializerSphere dummy(cx, cy, cz, radius, value);
		putput.addComponent(dummy.asString());
	}
	else{
		return 0;
	}

	// create diffusion tensor
	//std::shared_ptr<TensorConstant> DConstant = std::make_shared<TensorConstant>(rank, Dw);
	std::shared_ptr<TensorSpatial2D> DSpatial2D;
	std::shared_ptr<TensorSpatial3D> DSpatial3D;
	std::shared_ptr<dolfin::Expression> D;
	if(dimensions == 2){
		DSpatial2D = std::make_shared<TensorSpatial2D>(rank, Dw, Dg, get_10on10_test_cm());
		putput.addComponent(DSpatial2D->asString());
		D = DSpatial2D;
	}
	else if(dimensions == 3){
		DSpatial3D = std::make_shared<TensorSpatial3D>(rank, Dw, Dg, get_10on10on10_test_cm());
		putput.addComponent(DSpatial3D->asString());
		D = DSpatial3D;
	}
	else{
		return 0;
	}


	// create pde problem
	std::shared_ptr<ReactionDiffusionProblem> problem = std::make_shared<ReactionDiffusionProblem>(rank, mesh, D, rho, dt, theta);
	putput.addComponent(problem->asString());

	// create solver
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>();
	solver->parameters["linear_solver"] = "lu";
	solver->parameters["convergence_criterion"] = "incremental";
	solver->parameters["maximum_iterations"] = 50;
	solver->parameters["relative_tolerance"] = 1e-10;
	solver->parameters["absolute_tolerance"] = 1e-10;

	// create time stepper
	TimeStepper timeStepper = TimeStepper(rank, problem, solver, T, dt_min, dt_max, richTol, richSafe);
	putput.addComponent(timeStepper.asString());

	// create output files
	// simulation output
	auto outPvdReturn = putput.getFilePath(tagFolder, tagFile, "pvd");
	if(!outPvdReturn.first){	// check if path to file loaded successful
		return 0;
	}
	std::shared_ptr<dolfin::File> file = std::make_shared<dolfin::File>(outPvdReturn.second);
	std::stringstream ss;
	ss << tagCsv << "-" << rank;	// rank specific csv file
	auto outCsvReturn = putput.getFilePath(tagFolder, ss.str(), "csv");
	if(!outCsvReturn.first){	// check if path to file loaded successful
		return 0;
	}

	// create pre-simulation info
	putput.createRunInfo(tagFolder, tagFile);

	// run simulation
	RuntimeTracker tracker3 = timeStepper.run(timeAdaption, verbose, initialCondition,
			file, outCsvReturn.second,
			framesPerTimeUnit, dt);

	// overwrite INFO file with post-simulation details
	putput.addComponent(tracker3.asString());
	putput.createRunInfo(tagFolder, tagFile);

	MPI_Finalize(); //seems to trigger an abort
}



