
#include <vector>
#include <dolfin.h>
#include <Eigen/Dense>
#include <fstream>

#include "ProblemSolverContainer.h"
#include "FisherNewtonContainer.h"
#include "TimeStepper.h"
#include "RuntimeTracker.h"
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
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

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
	double dt_init = atof(argv[13]);
	double dt_max = atof(argv[14]);
	double T = atof(argv[15]);
	int framesPerTimeUnit = atoi(argv[16]);
	// reaction-diffusion coefficient
	////////////////////////////////////
	double Dw = atof(argv[17]);
	double Dg = atof(argv[18]);
	double rho = atof(argv[19]);
	double theta = atof(argv[20]);
	// timestepping
	///////////////////////////////////
	int verbose = atoi(argv[21]);
	int timeAdaption = atoi(argv[22]);
	double richTol = atof(argv[23]);
	double richSafe = atof(argv[24]);
	// solver
	///////////////////////////////////
	double newtontolrel = atof(argv[25]);
	double newtontolabs = atof(argv[26]);
	int newtonmaxiter = 50;
	double krylovtolrel = atof(argv[27]);
	double krylovtolabs = atof(argv[28]);
	int krylovmaxiter = 50;
	std::string ls = argv[29];
	std::string pc = argv[30];

	if(rank == 0 && verbose > 3){ std::cout << "load components..." << std::endl; };

	// in/output handler
	ReaderWriter putput = ReaderWriter(rank, outputParent, meshParent);

	// mesh read-in
	if(rank == 0 && verbose > 3){ std::cout << "	read in mesh..." << std::endl; };
	std::shared_ptr<dolfin::Mesh> mesh;
	std::pair<std::string, std::string> meshInfo = putput.loadMesh(meshName);
	if(meshInfo.first == "h5"){
		mesh = std::make_shared<dolfin::Mesh>();
		auto hdf5 = dolfin::HDF5File(MPI_COMM_WORLD, meshInfo.second, std::string("r"));
		hdf5.read(*mesh, "/mesh", false);
		// add mesh to components list before returning
		std::stringstream ss;
		ss << "Mesh: " << meshName << " (" << mesh->geometry().dim() << "D)" << std::endl;
		putput.addComponent(ss.str());
	}
	else if(meshInfo.first == "xml"){
		mesh = std::make_shared<dolfin::Mesh>(meshInfo.second);
		// add mesh to components list before returning
		std::stringstream ss;
		ss << "Mesh: " << meshName << " (" << mesh->geometry().dim() << "D)" << std::endl;
		putput.addComponent(ss.str());
	}
	else{
		if(rank == 0){
			std::cout << "WARNING (main) failed reading in mesh..." << std::endl;
			return 0;
		}
	}
	int dimensions = mesh->geometry().dim();
	if(rank == 0 && verbose > 3){ std::cout << "	mesh loaded!" << std::endl; };

	// create initial condition
	std::shared_ptr<dolfin::Expression> initialCondition;
	if(dimensions == 2){ initialCondition = std::make_shared<InitializerCircle>(cx, cy, radius, value);	}
	else if(dimensions == 3){ initialCondition = std::make_shared<InitializerSphere>(cx, cy, cz, radius, value); }
	else{ return 0; }
	std::stringstream ss;
	ss << "InitialCondition: " << (dimensions == 2 ? "Circle" : "Sphere") << " at (" << cx << ", " << cy << ", " << cz << "), r= " << radius << ", v= " << value << std::endl;
	putput.addComponent(ss.str());
	if(rank == 0 && verbose > 3){ std::cout << "	initial condition loaded!" << std::endl; };

	// create diffusion tensor
	std::shared_ptr<dolfin::Expression> D;
	if(dimensions == 2){ D = std::make_shared<TensorSpatial2D>(rank, Dw, Dg, get_10on10_test_cm());	}
	else if(dimensions == 3){ D = std::make_shared<TensorSpatial3D>(rank, Dw, Dg, get_10on10on10_test_cm()); }
	else{ return 0;	}
	std::stringstream ss2;
	ss2 << "Diffusion Tensor: D_white= " << Dw << ", D_grey= " << Dg << std::endl;
	putput.addComponent(ss2.str());
	if(rank == 0 && verbose > 3){ std::cout << "	D tensor loaded!" << std::endl; };


	// create FisherNewtonContainter
	FisherNewtonContainer problemContainer = FisherNewtonContainer(rank,
			mesh, initialCondition, D, rho, theta, dt_init);
	problemContainer.initializeSolver((verbose > 2 ? 1 : 0), newtontolrel, newtontolabs, newtonmaxiter,
				krylovtolrel, krylovtolabs, krylovmaxiter, ls, pc);
	std::stringstream ss3;
	ss3 << "Containter:	Mesh, Initial condition, Diffusion Tensor from above" << std::endl <<
			"		Reaction coefficient= " << rho << ", theta= " << theta << ", dt_init= " << dt_init << std::endl;
	putput.addComponent(ss3.str());
	if(rank == 0 && verbose > 3){ std::cout << "	ProblemContainer loaded!" << std::endl; };

	// create time stepper
	TimeStepper timeStepper = TimeStepper(rank, dt_min, dt_max, richTol, richSafe);
	putput.addComponent(timeStepper.asString());
	if(rank == 0 && verbose > 3){ std::cout << "	Timestepper loaded!" << std::endl; };

	// create output files
	// simulation output
	auto outPvdReturn = putput.getFilePath(tagFolder, tagFile, "pvd");
	if(!outPvdReturn.first){ return 0; }
	std::shared_ptr<dolfin::File> pvdFile = std::make_shared<dolfin::File>(outPvdReturn.second);
	// iteration details output
	auto outCsvReturn = putput.getFilePath(tagFolder, tagCsv, "csv");
	if(!outCsvReturn.first){ return 0; }
	if(rank == 0 && verbose > 3){ std::cout << "	Output files loaded!" << std::endl; };
	std::stringstream ss4;
	ss4 << "Files: pvd= " << outPvdReturn.second << std::endl <<
			"	csv= " << outCsvReturn.second << std::endl;
	putput.addComponent(ss4.str());


	// create pre-simulation info, print added components
	if(rank == 0){ putput.createRunInfo(tagFolder, tagFile); }

	// run simulation
	RuntimeTracker tracker3 = timeStepper.run(timeAdaption, verbose,
			pvdFile, outCsvReturn.second, framesPerTimeUnit,
			&problemContainer, T, dt_init);

	// overwrite INFO file with post-simulation details
	putput.addComponent(tracker3.asString());
	putput.createRunInfo(tagFolder, tagFile);

	// output timings to text file
	std::ofstream timings;
	timings.open(putput.getFilePath(tagFolder, "timings-latex", "txt").second, std::ios_base::trunc);
	std::set<dolfin::TimingType> s = {dolfin::TimingType::wall};
	timings << dolfin::timings(dolfin::TimingClear::keep, s).str_latex();
	timings.close();

	MPI_Finalize(); //seems to trigger an abort
}



