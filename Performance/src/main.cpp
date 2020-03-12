
#include <dolfin.h>
#include <petscksp.h>
#include <petscsys.h>
#include <fstream>
#include "../../DolfinFisherSolver/Tensors.h"
#include "../../DolfinFisherSolver/Initializers.h"
#include "../../DolfinFisherSolver/FisherProblem.h"

std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> create_10on10_vm()
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
std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> create_10on10on10_vm(){
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vm_w;
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vm_g;

	for(int i = 0; i < 16; i++){
		auto ok = create_10on10_vm();
		vm_w.push_back(ok[0]);
		vm_g.push_back(ok[1]);
	}

	std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> vm;
	vm.push_back(vm_w);
	vm.push_back(vm_g);

	return vm;
}


void runNewton(int rank, int nprocs, std::shared_ptr<dolfin::Mesh> mesh, std::string ls, std::string pc, double tol, std::string restype){
	std::shared_ptr<dolfin::Expression> D3 = std::make_shared<TensorConstant>(rank, 0.013 ); //create_10on10on10_vm()
	std::shared_ptr<dolfin::Expression> init3D = std::make_shared<InitializerSphere>(0.5, 0.5, 0.5, 0.1, 1);
	std::shared_ptr<FisherProblem> problem = std::make_shared<FisherProblem>(rank, mesh, D3, 0.025, 0.0000001, 1);
	auto us = problem->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *init3D;
	*u = *init3D;

	std::stringstream iters;
	iters << "type-1-all-iterations-time-test-procs-" << nprocs << "-tol" << tol << "-" << restype << ".csv";
	std::ofstream iterfile(iters.str(), std::ios_base::app);

	if(rank==0){ std::cout << "Solve " << ls << " + " << pc << std::endl;}
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>();
	solver->parameters["error_on_nonconvergence"] = false;
	solver->parameters["convergence_criterion"] = restype;
	solver->parameters["linear_solver"] = ls;
	solver->parameters["preconditioner"] = pc;
	solver->parameters["relative_tolerance"] = tol;
	solver->parameters["absolute_tolerance"] = tol;
	solver->parameters("krylov_solver")["relative_tolerance"] = tol;
	solver->parameters("krylov_solver")["absolute_tolerance"] = tol;
	dolfin::Timer t("BBB solve " + ls + " + " + pc);
	auto r = solver->solve(*problem, *u->vector());
	t.stop();
	if(rank == 0){ iterfile << ls << " + " << pc << "," << r.first << "," << solver->krylov_iterations() << "," << std::get<0>(t.elapsed()) << "," <<
		solver->relative_residual() << "," << solver->residual() << std::endl;
		std::cout << ls << " + " << pc << "," << r.first << "," << solver->krylov_iterations() << "," << std::get<0>(t.elapsed()) << "," <<
				solver->relative_residual() << "," << solver->residual() << std::endl;
	}
	iterfile.close();
}

void runCG(int rank, int nprocs, std::shared_ptr<dolfin::Mesh> mesh, std::string pc, double tol, std::string restype){
	std::shared_ptr<dolfin::Expression> D3 = std::make_shared<TensorSpatial3D>(rank, 0.013, 0.0013, create_10on10on10_vm());
	std::shared_ptr<dolfin::Expression> init3D = std::make_shared<InitializerSphere>(0.5, 0.5, 0.5, 0.1, 1);
	std::shared_ptr<FisherProblem> problem = std::make_shared<FisherProblem>(rank, mesh, D3, 0.025, 0.0000001, 1);
	auto us = problem->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *init3D;
	*u = *init3D;

	std::stringstream iters;
	iters << "cg-convergence-" << nprocs << "-tol-" << tol << "-" << restype << ".csv";
	std::ofstream iterfile(iters.str(), std::ios_base::app);

	if(rank==0){ std::cout << "Solve cg with preconditioner " << pc << std::endl;}
	std::shared_ptr<dolfin::PETScKrylovSolver> krylov = std::make_shared<dolfin::PETScKrylovSolver>("cg", pc);
	//krylov->set_nonzero_guess(true);
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>(
			problem->getMesh()->mpi_comm(), krylov, dolfin::PETScFactory::instance());
	solver->set_relaxation_parameter(0.9);
	solver->parameters["error_on_nonconvergence"] = false;
	solver->parameters["convergence_criterion"] = restype;
	solver->parameters["relative_tolerance"] = tol;
	solver->parameters["absolute_tolerance"] = tol;
	solver->parameters("krylov_solver")["relative_tolerance"] = tol;
	solver->parameters("krylov_solver")["absolute_tolerance"] = tol;

	// setup krylov residual historys
	PetscInt size = 1000;
	PetscReal *a;
	PetscCalloc1(size, &a);
	KSPSetResidualHistory(krylov->ksp(), a, 1000, PETSC_TRUE);
	// solve
	dolfin::Timer t("BBB solve cg + " + pc);
	auto r = solver->solve(*problem, *u->vector());
	t.stop();
	// get residual data
	PetscInt used = 0;
	PetscReal *residuals;
	KSPGetResidualHistory(krylov->ksp(), &residuals, &used);
	if(rank == 0){ iterfile << "preconditioner, " << pc << std::endl <<
			"Newton rel residual, " << solver->relative_residual() << std::endl <<
			"Newton abs residual, " << solver->residual() << std::endl <<
			"Newton iterations," << r.first << std::endl <<
			"Krylov iterations, " << solver->krylov_iterations() << std::endl <<
			"Elapsed time, " << std::get<0>(t.elapsed()) << std::endl <<
			"Iteration residuals, ";
		for(int i = 0; i < used; i++){
			iterfile << residuals[i] << ",";
		}
		iterfile << std::endl << std::endl;
	}

	iterfile.close();
	PetscFree(a);
}

void runBug(int rank, int nprocs, std::shared_ptr<dolfin::Mesh> mesh, std::string ls, std::string pc, double tol, std::string restype, bool bug){
	std::stringstream iters;
	iters << "bug-report-krylov-tolerance" << ".txt";
	std::ofstream iterfile(iters.str(), std::ios_base::app);

	if(rank==0){
		std::cout << "Running krylov tolerance bug (" << (bug ? "on" : "off") << ") with " << ls << " - " << pc << " - " << restype << std::endl;
		iterfile << "krylov tolerance bug (" << (bug ? "on" : "off") << ") with " << ls << " - " << pc << " - " << restype << std::endl;
	}

	// setup problem
	std::shared_ptr<dolfin::Expression> D3 = std::make_shared<TensorConstant>(rank, 0.013);
	std::shared_ptr<dolfin::Expression> init3D = std::make_shared<InitializerSphere>(0.5, 0.5, 0.5, 0.1, 1);
	std::shared_ptr<FisherProblem> problem = std::make_shared<FisherProblem>(rank, mesh, D3, 0.025, 0.0000001, 1);
	auto us = problem->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *init3D;
	*u = *init3D;

	// instantiate solver
	std::shared_ptr<dolfin::PETScKrylovSolver> krylov = std::make_shared<dolfin::PETScKrylovSolver>(ls, pc);
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>(
			problem->getMesh()->mpi_comm(), krylov, dolfin::PETScFactory::instance());
	//KSPSetInitialGuessNonzero(krylov->ksp(),PETSC_TRUE);
	krylov->set_from_options();
	solver->parameters["error_on_nonconvergence"] = false;
	solver->parameters["convergence_criterion"] = restype;
	solver->parameters["relative_tolerance"] = tol;
	solver->parameters["absolute_tolerance"] = tol;

	// trigger bug by setting tolerance
	if(bug){
		solver->parameters("krylov_solver")["relative_tolerance"] = tol;
		solver->parameters("krylov_solver")["absolute_tolerance"] = tol;
	}

	// setup krylov residual historys
	PetscInt size = 1000;
	PetscReal *a;
	PetscCalloc1(size, &a);
	KSPSetResidualHistory(krylov->ksp(), a, 1000, PETSC_FALSE);

	// solve
	dolfin::Timer t("BBB solve experimental" + ls + " + " + pc + " + " + restype);
	auto r = solver->solve(*problem, *u->vector());
	t.stop();

	// get residual data
	PetscInt used = 0;
	PetscReal *residuals;
	KSPGetResidualHistory(krylov->ksp(), &residuals, &used);
	if(rank == 0){
		std::cout <<
			"Newton rel residual, " << solver->relative_residual() << std::endl <<
			"Newton abs residual, " << solver->residual() << std::endl <<
			"Newton iterations," << r.first << std::endl <<
			"Krylov iterations, " << solver->krylov_iterations() << std::endl <<
			"Elapsed time, " << std::get<0>(t.elapsed()) << std::endl <<
			"Iteration residuals, ";
		iterfile <<
			"Newton rel residual, " << solver->relative_residual() << std::endl <<
			"Newton abs residual, " << solver->residual() << std::endl <<
			"Newton iterations," << r.first << std::endl <<
			"Krylov iterations, " << solver->krylov_iterations() << std::endl <<
			"Elapsed time, " << std::get<0>(t.elapsed()) << std::endl <<
			"Iteration residuals, ";
		for(int i = 0; i < used; i++){
			std::cout << residuals[i] << ",";
			iterfile << residuals[i] << ",";
		}
		iterfile << std::endl << std::endl;
		std::cout << std::endl << std::endl;
	}

	iterfile.close();
	PetscFree(a);
}

// automatic performance assestment of FisherSolver
int main(int argc, char* argv[]){
	dolfin::SubSystemsManager::init_mpi();

	// Parse command line options (will intialise PETSc if any PETSc
	// options are present, e.g. --petsc.pc_type=jacobi)
	dolfin::parameters.parse(argc, argv);
	dolfin::SubSystemsManager::init_petsc();

	// Default parameters
	dolfin::Parameters application_parameters("application_parameters");
	application_parameters.add("scaling_type", "weak", {"weak", "strong"});
	application_parameters.add("ndofs", 500000);
	application_parameters.add("output", false);
	application_parameters.add("type", 2);
	application_parameters.add("dolflog", 30);
	application_parameters.add("tol", 0.00000001);
	application_parameters.add("residual_type", "incremental");
	application_parameters.add("ls", "cg");
	application_parameters.add("pc", "jacobi");


	// Update from command line
	application_parameters.parse(argc, argv);

	// extract parameters
	const std::string problem_type = "fisher";
	const std::string scaling_type = application_parameters["scaling_type"];
	const int ndofs = application_parameters["ndofs"];
	const bool output = application_parameters["output"];
	const std::string out_dir = "output/";
	const int type = application_parameters["type"];
	const int dolflog = application_parameters["dolflog"];
	const double tol = application_parameters["tol"];
	const std::string residual_type = application_parameters["residual_type"];
	const std::string ls = application_parameters["ls"];
	const std::string pc = application_parameters["pc"];

	 // Set mesh partitioner
	 dolfin::parameters["mesh_partitioner"] = "SCOTCH";

	 int nprocs = dolfin::MPI::size(MPI_COMM_WORLD);;
	 int rank = dolfin::MPI::rank(MPI_COMM_WORLD);

	// calculate dofs
	int totalDofs = ndofs;
	if(scaling_type == "weak"){ totalDofs = nprocs * ndofs; }

	// create and prepate 3D problem
	// 3D: tetrahedras and 1 quantity of interessed is assumed
	dolfin::Timer t1("AAA create mesh 3D");
	int nxyz = pow((1.0 * totalDofs / 1.30) / 6, 1.0/3.0);
	std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(nxyz, nxyz, nxyz));
	t1.stop();
	std::shared_ptr<dolfin::Expression> D3 = std::make_shared<TensorSpatial3D>(rank, 0.013, 0.0013, create_10on10on10_vm());
	std::shared_ptr<dolfin::Expression> init3D = std::make_shared<InitializerSphere>(0.5, 0.5, 0.5, 0.1, 1);
	std::shared_ptr<FisherProblem> problem = std::make_shared<FisherProblem>(rank, mesh, D3, 0.025, 0.0000001, 1);
	auto us = problem->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *init3D;
	*u = *init3D;

	// Print simulation summary to console
	if (rank == 0)
	{
		std::stringstream ss;
		ss << "----------------------------------------------------------------" << std::endl <<
		 "Test problem summary" << std::endl <<
		"  Problem type:   "   << (type == 1 ? "all" : (type == 2 ? "cg as linsolver" : "unknown")) << " (3D)" << std::endl <<
		"  Scaling type:   "   << scaling_type << std::endl <<
		"  Num processes:  "  << nprocs << std::endl <<
		"  Mesh elements:  " <<  mesh->num_cells() << " (3D)" << std::endl <<
		"  Total dof:      " << us.at(0)->function_space()->dim() << " (3D)" << std::endl <<
		"  Average dof per rank: " << us.at(0)->function_space()->dim()/dolfin::MPI::size(mesh->mpi_comm()) << " (3D)" << std::endl <<
		"----------------------------------------------------------------" << std::endl;
		std::cout << ss.str() << std::endl;
	}

	// different solver - preconditioner pairs
	///////////////////////////////////////////////////
	// NOTE:
	//		no parallel support for ilu/ icc preconditioner
	//		no support for petsc linear solver
	//		richardson and jacobi doesn't converge
	//
	dolfin::set_log_level(dolflog);

	// tpye 1: run all, track time and iterations
	if(type == 1){
		std::stringstream iters;
		iters << "type-1-all-iterations-time-test-procs-" << nprocs << "-tol" << tol << "-" << residual_type << ".csv";
		std::ofstream iterfile(iters.str(), std::ios_base::trunc);
		iterfile << "name, newton iterations, krylov iterations, time for solve(), newton relative residual, newton abs residual" << std::endl;
		iterfile.close();

		runNewton(rank, nprocs, mesh,  "gmres", "petsc_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "gmres", "jacobi", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "gmres", "hypre_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "gmres", "hypre_euclid", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "cg", "petsc_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "cg", "jacobi", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "cg", "hypre_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "cg", "hypre_euclid", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "richardson", "petsc_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "richardson", "hypre_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "richardson", "hypre_euclid", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "bicgstab", "hypre_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "bicgstab", "hypre_euclid", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "minres", "hypre_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "minres", "hypre_euclid", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "tfqmr", "hypre_amg", tol, residual_type);
		runNewton(rank, nprocs,  mesh,  "tfqmr", "hypre_euclid", tol, residual_type);
	}
	// type 2: test only preconditioners for cg
	else if(type == 2){
		std::stringstream iters;
		iters << "cg-convergence-" << nprocs << "-tol-" << tol << "-" << residual_type << ".csv";
		std::ofstream iterfile(iters.str(), std::ios_base::trunc);
		iterfile.close();

		runCG(rank, nprocs, mesh, "jacobi", tol, residual_type);
		runCG(rank, nprocs, mesh, "petsc_amg", tol, residual_type);
		runCG(rank, nprocs, mesh, "hypre_amg", tol, residual_type);
		runCG(rank, nprocs, mesh, "hypre_euclid", tol, residual_type);
	}
	// showing the bug when setting the tolerance of the PETSc krylov solver
	// can be run with any solver but "richardson", dofs don't matter
	// use --type 4
	else if(type == 4){
		std::stringstream iters;
		iters << "bug-report-krylov-tolerance" << ".txt";
		std::ofstream iterfile(iters.str(), std::ios_base::trunc);
		iterfile.close();

		runBug(rank, nprocs,  mesh,  ls, pc, tol, "incremental", true);
		runBug(rank, nprocs,  mesh,  ls, pc, tol, "incremental", false);
		runBug(rank, nprocs,  mesh,  ls, pc, tol, "residual", false);
		runBug(rank, nprocs,  mesh,  ls, pc, tol, "residual", true);
	}

	else{
		if(rank == 0){ std::cout << "don't know what type to run" << std::endl;}
	}

	 // Display timings
	 dolfin::list_timings(dolfin::TimingClear::keep, {dolfin::TimingType::wall});

	 // output timings to text file
	 if(output){
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
