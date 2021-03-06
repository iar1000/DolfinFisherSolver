
#include <dolfin.h>
#include <petscksp.h>
#include <petscsys.h>
#include <fstream>
#include <string>
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
		bool krylovZeroStart,
		bool useBuffer, bool preLoadBuffer, std::string bufferCsv,
		int quadrature_deg){

	// setup petsc logging
	PetscLogStage stage1, stage2, stage3, stage4;

	// initialize problem
	PetscLogStageRegister("Problem initialization", &stage1);
	PetscLogStagePush(stage1);
	std::shared_ptr<TensorSpatial3D> D = std::make_shared<TensorSpatial3D>(rank, diffCoef1, diffCoef2, brain->getConcentrationMap(-1).first, translation, useBuffer);
	std::shared_ptr<FisherProblem> problem = std::make_shared<FisherProblem>(rank, quadrature_deg, mesh, D, reactCoef, 0.0000001, 1);
	std::shared_ptr<dolfin::Expression> initial_condition = std::make_shared<InitializerSphere>(initial_coordinates.at(0), initial_coordinates.at(1), initial_coordinates.at(2), 3, 1);
	auto us = problem->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *initial_condition;
	*u = *initial_condition;

	// instantiate solver
	PetscLogStageRegister("Solver initialization", &stage2);
	PetscLogStagePush(stage2);
	std::shared_ptr<dolfin::PETScKrylovSolver> krylov = std::make_shared<dolfin::PETScKrylovSolver>(ls, pc);
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>(
			problem->getMesh()->mpi_comm(), krylov, dolfin::PETScFactory::instance());
	if(krylovZeroStart){ KSPSetInitialGuessNonzero(krylov->ksp(),PETSC_TRUE); }
	PetscLogStagePop();

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

	if(rank==0){ std::cout << "Solve " << ls << " + " << pc <<  " -" << diffCoef1 << "-" << diffCoef2 << "-" << reactCoef << std::endl;}

	// solve without buffer
	PetscLogStageRegister("Solve w/o buffer", &stage3);
	PetscLogStagePush(stage3);
	dolfin::Timer t1("AAA solve " + ls + " + " + pc + " w/o buffer");
	auto r = solver->solve(*problem, *u->vector());
	*u0->vector() = *u->vector();
	t1.stop();
	PetscLogStagePop();
	double krylov_iters_nobuff = solver->krylov_iterations();

	if(rank == 0) {
		std::cout << "	" <<  "w/o" << " buffer (" << (quadrature_deg == 6 ? "Q6" : "Q4") << "): " << std::get<0>(t1.elapsed()) << " (krylov : " << krylov_iters_nobuff << ", newton: " << r.first  << ") " << std::endl;
		std::cout << "	buffer miss count: " << D->getMissCount() << std::endl;
	}

	// dump everything to csv
	// D->dumpBuffer("bufferDump.csv");

	// solve with buffer
	PetscLogStageRegister("Solve w buffer", &stage4);
	PetscLogStagePush(stage4);
	dolfin::Timer t("AAA solve " + ls + " + " + pc + " w buffer");
	r = solver->solve(*problem, *u->vector());
	t.stop();
	PetscLogStagePop();

	if(rank == 0) {
			std::cout << "	w buffer (" << (quadrature_deg == 6 ? "Q6" : "Q4") << "): " << std::get<0>(t.elapsed()) << " (krylov : " << solver->krylov_iterations() << ", newton: " << r.first << ") " << std::endl;
			std::cout << "	buffer miss count: " << D->getMissCount() << std::endl;
	}

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


static std::vector<std::pair<std::string, std::string>> combis_cg{
	{"cg", "petsc_amg"}, {"cg", "hypre_amg"}, {"cg", "hypre_euclid"}, {"cg", "jacobi"}
};
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
    application_parameters.add("meshname", "lh-plial-dof-5500k.h5");
	application_parameters.add("type", 2);
	application_parameters.add("name", "no-name");
	application_parameters.add("newton_tol", 0.00000001);
	application_parameters.add("ls", "cg");
	application_parameters.add("pc", "hypre_euclid");
	application_parameters.add("diffCoef1", 0.013);
	application_parameters.add("diffCoef2", 0.0013);
	application_parameters.add("reactCoef", 0.025);
	application_parameters.add("krylovnonzero", false);
	application_parameters.add("useBuffer", 1);
	application_parameters.add("quadrature", 4);
	application_parameters.add("preloadBuffer", 0);				// not used
	application_parameters.add("bufferCsv", "bufferDump.csv");	// not used


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
	const int useBuffer = application_parameters["useBuffer"];
	const int quadratureDegree = application_parameters["quadrature"];
	const int preloadBuffer = application_parameters["preloadBuffer"];
	const std::string bufferCsv = application_parameters["bufferCsv"];


	// Set mesh partitioner
	 dolfin::parameters["mesh_partitioner"] = "SCOTCH";

	 int nprocs = dolfin::MPI::size(MPI_COMM_WORLD);;
	 int rank = dolfin::MPI::rank(MPI_COMM_WORLD);

	 // create mesh and dummy variables for intel
	 if(rank == 0){ std::cout << "setup performance test.." << std::endl; }
	 std::shared_ptr<dolfin::Mesh> mesh;
	 std::vector<std::string> tokens;
	 std::stringstream ss(meshname);
	 std::string token;
	 while(std::getline(ss, token, '.')){ tokens.push_back(token); }
	 if(tokens.at(1) == "h5"){
		 mesh = std::make_shared<dolfin::Mesh>();
		 auto hdf5 = dolfin::HDF5File(MPI_COMM_WORLD, mesh_path+meshname, std::string("r"));
		 hdf5.read(*mesh, "/mesh", false);
	 }
	 else if(tokens.at(1) == "xdmf"){
		 if(rank == 0){ std::cout << "\treading in xdmf mesh.." << std::endl; }
		 mesh = std::make_shared<dolfin::Mesh>(MPI_COMM_WORLD);
		 dolfin::XDMFFile(MPI_COMM_WORLD, mesh_path+meshname).read(*mesh);
	 }
	 else if(tokens.at(1) == "xml"){
		 mesh = std::make_shared<dolfin::Mesh>(mesh_path + meshname);
	 }
	 else { std::cout << "wrong mesh format!" << std::endl; return 0; }
	 std::shared_ptr<Brain> brain = std::make_shared<Brain>(rank, 2, "../../brain-data/brainweb");
	 std::shared_ptr<dolfin::Expression> dummy_D = std::make_shared<TensorConstant>(rank, 0.1);
	 std::shared_ptr<FisherProblem> dummy_problem = std::make_shared<FisherProblem>(rank, quadratureDegree, mesh, dummy_D, 0.1, 0.1, 0.1);
	 int ndofs = dummy_problem->getUs().at(0)->function_space()->dim();

	// Print simulation summary to console
	if (rank == 0)
	{
		std::stringstream ss;
		ss << "----------------------------------------------------------------" << std::endl <<
		 "Test problem summary" << std::endl <<
		"  Problem type:   "   << type << std::endl <<
		"  Num processes:  "  << nprocs << std::endl <<
		"  Mesh:           " << mesh_path+meshname << std::endl <<
		"  Mesh elements (approx.):  " <<  mesh->num_cells() * nprocs << std::endl <<
		"  Global dof:      " << ndofs << std::endl <<
		"  Average dof per rank: " << ndofs/nprocs << std::endl <<
		"----------------------------------------------------------------" << std::endl;
		std::cout << ss.str() << std::endl;
	}

	if(rank == 0){ std::cout << "run performance tests..." << std::endl; }
	dolfin::set_log_level(dolflog);

	// type 1: test all ls,pc combinations with standart parameters
	if(type == 1){
		std::stringstream iters;
		iters << name << "-type-1-tol-" << newton_tol << "-nprocs-" << nprocs << "-dofpr-" << (ndofs / nprocs) << ".csv";
		std::ofstream iterfile(out_dir + iters.str(), std::ios_base::app);
		if(rank == 0){ iterfile << output_format; }
		iterfile.close();

		// run test on all combinations
		for(int i = 0; i < combis.size(); i++){
			runTest(out_dir + iters.str(), rank, mesh, brain,
					combis.at(i).first, combis.at(i).second, newton_tol,
					diffCoef1, (diffCoef1/10), reactCoef, krylovnonzero, useBuffer ? true : false, preloadBuffer ? true : false, bufferCsv, quadratureDegree);
		}
	}
	// type 2: test conjugate gradient solver with hypre euclid preconditioner
	else if(type == 2){
		std::stringstream iters;
		iters << name << "-type-2-tol-" << newton_tol << "-nprocs-" << nprocs << "-dofpr-" << (ndofs / nprocs) << ".csv";
		std::ofstream iterfile(out_dir + iters.str(), std::ios_base::app);
		if(rank == 0){ iterfile << output_format; }
		iterfile.close();

		runTest(out_dir + iters.str(), rank, mesh, brain,
							"cg", "hypre_euclid", newton_tol,
							diffCoef1, (diffCoef1/10), reactCoef, krylovnonzero, useBuffer ? true : false, preloadBuffer ? true : false, bufferCsv, quadratureDegree);
	}
	// type 3: test conjugate gradient solver
	else if(type == 3){
		std::stringstream iters;
		iters << name << "-type-3-tol-" << newton_tol << "-nprocs-" << nprocs << "-dofpr-" << (ndofs / nprocs) << ".csv";
		std::ofstream iterfile(out_dir + iters.str(), std::ios_base::app);
		if(rank == 0){ iterfile << output_format; }
		iterfile.close();

		// run test on all combinations
		for(int i = 0; i < combis_cg.size(); i++){
			runTest(out_dir + iters.str(), rank, mesh, brain,
					combis_cg.at(i).first, combis_cg.at(i).second, newton_tol,
					diffCoef1, (diffCoef1/10), reactCoef, krylovnonzero, useBuffer ? true : false,preloadBuffer ? true : false, bufferCsv, quadratureDegree);
		}
	}
	// type 4: test bicgstab with hypre euclid preconditioner
	else if(type == 4){
		std::stringstream iters;
		iters << name << "-type-4-tol-" << newton_tol << "-nprocs-" << nprocs << "-dofpr-" << (ndofs / nprocs) << ".csv";
		std::ofstream iterfile(out_dir + iters.str(), std::ios_base::app);
		if(rank == 0){ iterfile << output_format; }
		iterfile.close();

		runTest(out_dir + iters.str(), rank, mesh, brain,
				"bicgstab", "hypre_euclid", newton_tol,
				diffCoef1, (diffCoef1/10), reactCoef, krylovnonzero, useBuffer ? true : false, preloadBuffer ? true : false, bufferCsv, quadratureDegree);
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

	return 0;
}
