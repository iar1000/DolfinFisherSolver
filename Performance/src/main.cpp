
#include <dolfin.h>
#include <petscksp.h>
#include <petscsys.h>
#include <fstream>
#include "../../DolfinFisherSolver/Tensors.h"
#include "../../DolfinFisherSolver/Initializers.h"
#include "../../DolfinFisherSolver/FisherProblem.h"
#include "../../DolfinFisherSolver/Brain.h"
#include "../../DolfinFisherSolver/FisherNewtonContainer.h"


// static variables
static std::string output_format = "ls, pc, dw, dg, rho, newton iterations, krylov iterations, time for solve(), newton relative residual, newton abs residual, residuals \n";
static const std::vector<int> translation = {20, 23, 33};
static const std::vector<int> initial_coordinates = {20, 90, 65};
static const std::string out_dir = "../output/";
static const std::string mesh_path = "../../mesh/";


void runTest(std::string filepath, int rank,
		std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<Brain> brain,
		std::string ls, std::string pc, double newton_tol,
		double diffCoef1, double diffCoef2, double reactCoef,
		bool krylovZeroStart){
	// initialize problem
	std::shared_ptr<dolfin::Expression> D = std::make_shared<TensorSpatial3D>(rank, diffCoef1, diffCoef2, brain->getConcentrationMap().first, translation);
	//std::shared_ptr<dolfin::Expression> D = std::make_shared<TensorConstant>(rank, diffCoef1);
	std::shared_ptr<FisherProblem> problem = std::make_shared<FisherProblem>(rank, mesh, D, reactCoef, 0.0000001, 1);
	std::shared_ptr<dolfin::Expression> initial_condition = std::make_shared<InitializerSphere>(initial_coordinates.at(0), initial_coordinates.at(1), initial_coordinates.at(2), 3, 1);
	auto us = problem->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *initial_condition;
	*u = *initial_condition;

	if(rank==0){ std::cout << "Solve " << ls << " + " << pc <<  " -" << diffCoef1 << "-" << diffCoef2 << "-" << reactCoef << std::endl;}

	// instantiate solver
	std::shared_ptr<dolfin::PETScKrylovSolver> krylov = std::make_shared<dolfin::PETScKrylovSolver>(ls, pc);
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>(
			problem->getMesh()->mpi_comm(), krylov, dolfin::PETScFactory::instance());
	if(krylovZeroStart){ KSPSetInitialGuessNonzero(krylov->ksp(),PETSC_TRUE); }

	// set solver parameters
	krylov->set_from_options();
	solver->parameters["error_on_nonconvergence"] = false;
	solver->parameters["convergence_criterion"] = "residual";
	solver->parameters["relative_tolerance"] = newton_tol/100;
	solver->parameters["absolute_tolerance"] = newton_tol;

	// setup krylov residual historys
	PetscInt size = 1000;
	PetscReal *a;
	PetscCalloc1(size, &a);
	KSPSetResidualHistory(krylov->ksp(), a, 1000, PETSC_FALSE);

	// solve
	dolfin::Timer t("AAA solve " + ls + " + " + pc);
	auto r = solver->solve(*problem, *u->vector());
	t.stop();

	// get residual data
	PetscInt used = 0;
	PetscReal *residuals;
	KSPGetResidualHistory(krylov->ksp(), &residuals, &used);


	// output to file
	if(rank == 0){
		// output file
		std::ofstream iterfile(filepath, std::ios_base::app);

		iterfile << ls << "," << pc <<  "," <<
		diffCoef1 << "," << diffCoef2 << "," << reactCoef << "," <<
		r.first << "," <<
		solver->krylov_iterations() << "," <<
		std::get<0>(t.elapsed()) << "," <<
		solver->relative_residual() << "," <<
		solver->residual();
		for(int i = 0; i < used; i++){
				iterfile << "," << residuals[i];
		}
		iterfile << std::endl;
		iterfile.close();
	}
}

void runConvergenceStudy(int rank, double dt){
	std::shared_ptr<Brain> brain = std::make_shared<Brain>(rank, 2, "../../brain-data/brainweb");
	std::shared_ptr<dolfin::Expression> initial_condition = std::make_shared<InitializerSphere>(initial_coordinates.at(0), initial_coordinates.at(1), initial_coordinates.at(2), 3, 1);
	std::shared_ptr<dolfin::Expression> D = std::make_shared<TensorSpatial3D>(rank, 0.13, 0.013, brain->getConcentrationMap().first, translation);

	//std::shared_ptr<dolfin::Mesh> mesh80 = std::make_shared<dolfin::Mesh>("../../mesh/lh-white-hull-flood-0-1-merge-5-dof-80k.xml");
	//std::shared_ptr<dolfin::Mesh> mesh600 = std::make_shared<dolfin::Mesh>("../../mesh/lh-white-hull-flood-0-1-merge-5-dof-600k.xml");
	//std::shared_ptr<dolfin::Mesh> mesh4700 = std::make_shared<dolfin::Mesh>("../../mesh/lh-white-hull-flood-0-1-merge-5-dof-4700k.xml");
	//std::shared_ptr<dolfin::Mesh> mesh37000 = std::make_shared<dolfin::Mesh>("../../mesh/lh-white-hull-flood-0-1-merge-5-dof-37000k.xml");

	int numMesh = 2;
	std::shared_ptr<dolfin::Mesh> mesh80 = std::make_shared<dolfin::Mesh>("../../mesh/box-10on10on10-res-23.xml");
	std::shared_ptr<dolfin::Mesh> mesh600 = std::make_shared<dolfin::Mesh>("../../mesh/box-10on10on10-res-29.xml");
	//std::shared_ptr<dolfin::Mesh> mesh4700 = std::make_shared<dolfin::Mesh>("../../mesh/box-10on10on10-res-40.xml");
	//std::shared_ptr<dolfin::Mesh> mesh37000 = std::make_shared<dolfin::Mesh>("../../mesh/box-10on10on10-res-50.xml");

	// create fisher containers
	std::vector<FisherNewtonContainer> containers;
	containers.push_back(FisherNewtonContainer(rank,	mesh80, initial_condition, D, 0.025, 1, dt));
	containers.push_back(FisherNewtonContainer(rank,	mesh600, initial_condition, D, 0.025, 1, dt));
			//FisherNewtonContainer(rank,	mesh4700, initial_condition, D, 0.025, 1, dt),
			//FisherNewtonContainer(rank,	mesh37000, initial_condition, D, 0.025, 1, dt)
	for(int i = 0; i < numMesh; i++){
		containers.at(i).initializeSolver(0, 0.00000001, 0.0000000001, 50,
					0, 0, 50, "cg", "jacobi");
	}

	// get coordinates of lowest resolution mesh
	std::vector<double> coordVec = mesh80->coordinates();
	double* coordsArr = &coordVec[0];
	int numPoints = coordVec.size() / 3;
	dolfin::Array<double> coordinates(coordVec.size(), coordsArr);
	dolfin::Array<double> values(numPoints);
	std::cout << rank << " has #coordinates= " << coordVec.size() << std::endl;
	double T = 10;
	double t = 0;

	while(t < T){
		if(rank == 0){ std::cout << "run time t= " << t << std::endl; }
		for(int i = 0; i < numMesh; i++){
			if(rank == 0){ std::cout << "	solve mesh size  " << i << std::endl; }

			// run iteration, continue on successful convergence
			int converged = containers.at(i).solve(t, dt);
			if(!converged){
				if(rank == 0){ std::cout << std::endl << "iteration failed to converge at t = " << t << std::endl; }
			}
		}
		// update t
		t += dt;

		if(rank == 0){ std::cout << "evaluate t= " << t << std::endl; }
		for(int i = 0; i < numMesh; i++){
			containers.at(i).getProblem()->getUs().at(0)->eval(values, coordinates);
			std::ofstream csv;
			std::stringstream name;
			name << "mesh-" << i << "-t-" << t << "evaluations.csv";
			csv.open(name.str(), std::ios_base::app);
			for(int j = 0; j < coordVec.size(); j += 3){
				csv << coordsArr[j] << "," << coordsArr[j+1] << "," << coordsArr[j+2] << "," << values[j/3] << std::endl;
			}
			csv.close();
		}

	}
}


static std::vector<std::pair<std::string, std::string>> combis_exc_euclid{
	{"gmres", "petsc_amg"}, {"gmres", "hypre_amg"}, {"gmres", "jacobi"},
	{"cg", "petsc_amg"}, {"cg", "hypre_amg"}, {"cg", "jacobi"},
	{"richardson", "petsc_amg"}, {"richardson", "hypre_amg"},
	{"bicgstab", "hypre_amg"},
	{"minres", "hypre_amg"},
	{"tfqmr", "hypre_amg"}
};
static std::vector<std::pair<std::string, std::string>> combis{
	{"gmres", "petsc_amg"}, {"gmres", "hypre_amg"}, {"gmres", "hypre_euclid"}, {"gmres", "jacobi"},
	{"cg", "petsc_amg"}, {"cg", "hypre_amg"}, {"cg", "hypre_euclid"}, {"cg", "jacobi"},
	{"richardson", "petsc_amg"}, {"richardson", "hypre_amg"}, {"richardson", "hypre_euclid"},
	{"bicgstab", "hypre_amg"}, {"bicgstab", "hypre_euclid"},
	{"minres", "hypre_amg"}, {"minres", "hypre_euclid"},
	{"tfqmr", "hypre_amg"}, {"tfqmr", "hypre_euclid"}
};

// automatic performance assestment of FisherSolver
int main(int argc, char* argv[]){
	dolfin::SubSystemsManager::init_mpi();

	// Parse command line options (will intialise PETSc if any PETSc
	// options are present, e.g. --petsc.pc_type=jacobi)
	dolfin::parameters.parse(argc, argv);
	dolfin::SubSystemsManager::init_petsc();

	// Default parameters
	dolfin::Parameters application_parameters("application_parameters");
	application_parameters.add("dolflog", 30);
	application_parameters.add("meshname", "lh-white-hull-flood-0-1-merge-5-dof-unknown.xml");
	application_parameters.add("type", 1);
	application_parameters.add("name", "no-name");
	application_parameters.add("newton_tol", 0.00000001);
	application_parameters.add("ls", "cg");
	application_parameters.add("pc", "jacobi");
	application_parameters.add("diffCoef1", 0.013);
	application_parameters.add("diffCoef2", 0.0013);
	application_parameters.add("reactCoef", 0.025);
	application_parameters.add("krylovnonzero", false);


	// Update from command line
	application_parameters.parse(argc, argv);

	// extract parameters
	const std::string meshname = application_parameters["meshname"];
	const int type = application_parameters["type"];
	const std::string name = application_parameters["name"];
	const int dolflog = application_parameters["dolflog"];
	const double newton_tol = application_parameters["newton_tol"];
	const std::string ls = application_parameters["ls"];
	const std::string pc = application_parameters["pc"];
	const double diffCoef1 = application_parameters["diffCoef1"];
	const double diffCoef2 = application_parameters["diffCoef2"];
	const double reactCoef = application_parameters["reactCoef"];
	const bool krylovnonzero = application_parameters["krylovnonzero"];

	 // Set mesh partitioner
	 dolfin::parameters["mesh_partitioner"] = "SCOTCH";

	 int nprocs = dolfin::MPI::size(MPI_COMM_WORLD);;
	 int rank = dolfin::MPI::rank(MPI_COMM_WORLD);

	 // create mesh and dummy variables for intel
	 std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>(mesh_path + meshname);
	 std::shared_ptr<Brain> brain = std::make_shared<Brain>(rank, 2, "../../brain-data/brainweb");
	 std::shared_ptr<dolfin::Expression> dummy_D = std::make_shared<TensorConstant>(rank, 0.1);
	 std::shared_ptr<FisherProblem> dummy_problem = std::make_shared<FisherProblem>(rank, mesh, dummy_D, 0.1, 0.1, 0.1);
	 int ndofs = dummy_problem->getUs().at(0)->function_space()->dim();

	// Print simulation summary to console
	if (rank == 0)
	{
		std::stringstream ss;
		ss << "----------------------------------------------------------------" << std::endl <<
		 "Test problem summary" << std::endl <<
		"  Problem type:   "   << type << std::endl <<
		"  Num processes:  "  << nprocs << std::endl <<
		"  Mesh:           " << mesh_path << std::endl <<
		"  Mesh elements:  " <<  mesh->num_cells() << std::endl <<
		"  Global dof:      " << ndofs << std::endl <<
		"  Average dof per rank: " << ndofs/nprocs << std::endl <<
		"----------------------------------------------------------------" << std::endl;
		std::cout << ss.str() << std::endl;
	}

	dolfin::set_log_level(dolflog);

	// type 1: test all ls,pc combinations with standart parameters
	if(type == 1){
		std::stringstream iters;
		iters << name << "-type-1-tol-" << newton_tol << "-nprocs-" << nprocs << "-dofpr-" << (ndofs / nprocs) << ".csv";
		std::ofstream iterfile(out_dir + iters.str(), std::ios_base::app);
		iterfile << output_format;
		iterfile.close();

		// run test on all combinations
		for(int i = 0; i < combis.size(); i++){
			runTest(out_dir + iters.str(), rank, mesh, brain,
					combis.at(i).first, combis.at(i).second, newton_tol,
					diffCoef1, (diffCoef1/10), reactCoef, krylovnonzero);
		}
	}
	else if(type == 3){
		runConvergenceStudy(rank, 0.00001);
	}

	else{
		if(rank == 0){ std::cout << "don't know what type to run" << std::endl;}
	}

 // output timings to text file
	 if(false){
		 std::ofstream timings;
		 std::stringstream ss;
		 ss << "timings-" << nprocs << ".txt";
		 timings.open(ss.str(), std::ios_base::trunc);
		 std::set<dolfin::TimingType> s = {dolfin::TimingType::wall};
		 timings << dolfin::timings(dolfin::TimingClear::keep, s).str_latex();
		 timings.close();

		 std::stringstream sss;
		 sss << "timings-" << nprocs << ".xml";
		 dolfin::dump_timings_to_xml(sss.str(), dolfin::TimingClear::keep);
	 }

	MPI_Finalize();
}
