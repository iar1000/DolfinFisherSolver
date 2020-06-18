
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
#include "Brain.h"

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

	for(int i = 0; i < 11; i++){
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
	int dt_runlength = atoi(argv[25]);		// passed into Timestepper.run(), but not used
	// solver
	///////////////////////////////////
	double newtontolrel = atof(argv[26]);
	double newtontolabs = atof(argv[27]);
	int newtonmaxiter = 50;
	double krylovtolrel = atof(argv[28]);	// passed into Problem-Solver container, but not used
	double krylovtolabs = atof(argv[29]);	// passed into Problem-Solver container, but not used
	int krylovmaxiter = 50;					// passed into Problem-Solver container, but not used
	std::string ls = argv[30];
	std::string pc = argv[31];
	// translation of mesh values to brainweb concentration values
	/////////////////////////////////////
	int calcTrans = atoi(argv[32]);
	std::vector<int> translation = {atoi(argv[33]), atoi(argv[34]), atoi(argv[35])};
	int sliceNumber = 60;
	std::string brainweb_path = argv[36];

	/* DEBUG INPUT OF BASH SCRIPTS
	for(int i = 0; i < 36; i++){
		std::cout << "argument " << i << " : " << argv[i] << std::endl;
	}
	return 0;
	*/

	if(rank == 0 && verbose > 3){ std::cout << "load components..." << std::endl; };

	// in/output handler
	ReaderWriter putput = ReaderWriter(rank, outputParent, meshParent);

	// status: simulation started
	if(rank == 0){ putput.updateStatusFile(0, tagFolder, nprocs); }

	// status: before mesh loading
	if(rank == 0){ putput.updateStatusFile(1, tagFolder, nprocs); }

	// mesh read-in
	if(rank == 0 && verbose > 3){ std::cout << "	read in mesh..." << std::endl; };
	std::shared_ptr<dolfin::Mesh> mesh;
	std::pair<std::string, std::string> meshInfo = putput.loadMesh(meshName);
	if(meshInfo.first == "xdmf"){
		mesh = std::make_shared<dolfin::Mesh>(MPI_COMM_WORLD);
		dolfin::XDMFFile(MPI_COMM_WORLD, meshInfo.second).read(*mesh);
	}
	else if(meshInfo.first == "h5"){
		mesh = std::make_shared<dolfin::Mesh>(MPI_COMM_WORLD);
		auto hdf5 = dolfin::HDF5File(MPI_COMM_WORLD, meshInfo.second, std::string("r"));
		hdf5.read(*mesh, "/mesh", false);
	}
	else if(meshInfo.first == "xml"){
		mesh = std::make_shared<dolfin::Mesh>(meshInfo.second);
	}
	else{
		if(rank == 0){
			std::cout << "WARNING (main) failed reading in mesh..." << std::endl;
			return 0;
		}
	}

	// add mesh to components list before returning
	std::stringstream meshstring;
	meshstring << "Mesh: " << meshName << " (" << mesh->geometry().dim() << "D)" << std::endl;
	putput.addComponent(meshstring.str());
	int dimensions = mesh->geometry().dim();
	if(rank == 0 && verbose > 3){ std::cout << "	mesh (" << mesh->geometry().dim() << "D) loaded!" << std::endl; };


	// create brain
	if(rank == 0 && verbose > 3){ std::cout << "	read in brain from " << brainweb_path << std::endl; };
	Brain brain(rank, verbose, brainweb_path);
	// calculate translation and return
	if(dimensions == 3 && calcTrans){
		if(nprocs > 1 && rank == 0){
			std::cout << "To calculate correct translation data, please run with only 1 processor. currently " << nprocs <<
					"\n\nSimulation stopped, change number processors to 1 to proceed";
			return 0;
		}
		auto optTransReturn = brain.greedyOptimalTranslation(mesh, 10, true, splitString(meshName, '.').at(0));
		translation = {optTransReturn.at(2), optTransReturn.at(3), optTransReturn.at(4)};

		std::stringstream ss4;
		ss4 << "Translation Coordinates= " << translation.at(0) << ", " << translation.at(1) << ", " << translation.at(2) << std::endl <<
					"Mesh to brainweb miss rate in %= " << (int)((double)optTransReturn.at(1) / (double)optTransReturn.at(0) * 100) << std::endl <<
					"\n\nRun intended to calculate translation data, no Fisher simulation run!\n\n";
		putput.addComponent(ss4.str());
		std::cout << ss4.str() << std::endl;
		// create pre-simulation info, print added components
		if(rank == 0){ putput.createRunInfo(tagFolder, tagFile); }

		return 0;
	}
	// get concentration map
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*> brainwebConcentrationMap;
	if(dimensions == 3){ brainwebConcentrationMap = brain.getConcentrationMap(-1); }
	else if(dimensions == 2){
		// check if desired slice is in bounds
		if(sliceNumber >= 0 && sliceNumber < 181){
			brainwebConcentrationMap = brain.getConcentrationMap(sliceNumber);
			// no translation necessary since mesh created directly from brainweb
			for(int i = 0; i < 3; i++){
				translation.at(i) = 0;
			}
		}
		else{ if(rank == 0){ std::cout << "desired slice is not in bounds of brainweb atlas!" << std::endl; }}
	}
	putput.addComponent(brain.asString());

	// create diffusion tensor
	std::shared_ptr<dolfin::Expression> D;
	if(dimensions == 2){ D = std::make_shared<TensorSpatial2D>(rank, Dw, Dg, get_10on10_test_cm());	}			// @ No value map yet
	else if(dimensions == 3){ D = std::make_shared<TensorSpatial3D>(rank, Dw, Dg, brainwebConcentrationMap.first, translation, true); }
	else{ return 0;	}
	std::stringstream ss2;
	ss2 << "Diffusion Tensor: D_white= " << Dw << ", D_grey= " << Dg << std::endl
			<< "Reaction Coefficient= " << rho << std::endl
			<< "Translation= " << translation.at(0) << ", " << translation.at(1) << ", " << translation.at(2) << std::endl;
	putput.addComponent(ss2.str());
	if(rank == 0 && verbose > 3){ std::cout << "	D tensor loaded!" << std::endl; };

	// create initial condition
	std::shared_ptr<dolfin::Expression> initialCondition;
	if(dimensions == 2){ initialCondition = std::make_shared<InitializerCircle>(cx, cy, radius, value);	}
	else if(dimensions == 3){ initialCondition = std::make_shared<InitializerGaussianSphere>(cx, cy, cz, radius, value); }
	else{ return 0; }
	std::stringstream ss;
	ss << "InitialCondition: " << (dimensions == 2 ? "Circle" : "GaussianSphere") << " at (" << cx << ", " << cy << ", " << cz << "), r= " << radius << ", v= " << value << std::endl;
	putput.addComponent(ss.str());
	if(rank == 0 && verbose > 3){ std::cout << "	initial condition loaded!" << std::endl; };


	// create FisherNewtonContainter
	FisherNewtonContainer problemContainer = FisherNewtonContainer(rank,
			mesh, initialCondition, D, rho, theta, dt_init); // first D = initialCondition
	problemContainer.initializeSolver((verbose > 3 ? 1 : 0), newtontolrel, newtontolabs, newtonmaxiter,
				krylovtolrel, krylovtolabs, krylovmaxiter, ls, pc);
	putput.addComponent(problemContainer.asString());
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
	if(rank == 0){
		putput.createRunInfo(tagFolder, tagFile);
		putput.updateStatusFile(2, tagFolder, nprocs);
	}

	// run simulation
	RuntimeTracker tracker3 = timeStepper.run(timeAdaption, verbose,
			pvdFile, outCsvReturn.second, framesPerTimeUnit,
			&problemContainer, T, dt_init, dt_runlength);

	// overwrite INFO file with post-simulation details
	if(rank == 0){
		putput.addComponent(tracker3.asString());
		putput.createRunInfo(tagFolder, tagFile);
		putput.updateStatusFile(3, tagFolder, nprocs);
	}

	// output timings to text file
	std::ofstream timings;
	timings.open(putput.getFilePath(tagFolder, "timings-latex", "txt").second, std::ios_base::trunc);
	std::set<dolfin::TimingType> s = {dolfin::TimingType::wall};
	timings << dolfin::timings(dolfin::TimingClear::keep, s).str_latex();
	timings.close();

	return 0;
}



