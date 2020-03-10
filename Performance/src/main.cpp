
#include <dolfin.h>
#include <fstream>
#include "../../DolfinFisherSolver/Tensors.h"
#include "../../DolfinFisherSolver/Initializers.h"
#include "../../DolfinFisherSolver/FisherProblem.h"
#include "../../DolfinFisherSolver/FisherNewtonContainer.h"

void runNewton(int rank, int nprocs, std::shared_ptr<dolfin::Function> u, std::shared_ptr<FisherProblem> problem,
		std::string ls, std::string pc){
	std::stringstream iters;
	iters << "solver-iterations-" << nprocs << ".txt";
	std::ofstream iterfile(iters.str(), std::ios_base::app);

	if(rank==0){ std::cout << "Solve " << ls << " + " << pc << std::endl;}
	dolfin::Timer t("BBB solve " + ls + " + " + pc);
	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>();
	solver->parameters["error_on_nonconvergence"] = false;
	solver->parameters["convergence_criterion"] = "incremental";
	solver->parameters["linear_solver"] = ls;
	solver->parameters["preconditioner"] = pc;
	auto r = solver->solve(*problem, *u->vector());
	t.stop();
	iterfile << ls << " + " << pc << std::endl <<
			"Newton iterations: " << r.first << std::endl <<
			"Krylov iterations: " << solver->krylov_iterations() << std::endl <<
			"Elapsed time: " << std::get<0>(t.elapsed()) << std::endl << std::endl;

	iterfile.close();
}

// automatic performance assestment of FisherSolver
int main(int argc, char* argv[]){
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	/* CHeck available solvers
	if(rank==0){
		dolfin::list_krylov_solver_preconditioners();
		dolfin::list_linear_solver_methods();
	}
	*/

	// Parse command line options (will intialise PETSc if any PETSc
	// options are present, e.g. --petsc.pc_type=jacobi)
	dolfin::parameters.parse(argc, argv);
	dolfin::SubSystemsManager::init_petsc();

	// Default parameters
	dolfin::Parameters application_parameters("application_parameters");
	application_parameters.add("scaling_type", "weak", {"weak", "strong"});
	application_parameters.add("ndofs", 500000);
	application_parameters.add("output", false);
	application_parameters.add("type", 1);
	application_parameters.add("dolflog", 30);

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

	 // Set mesh partitioner
	 dolfin::parameters["mesh_partitioner"] = "SCOTCH";

	// calculate dofs
	int totalDofs = ndofs;
	if(scaling_type == "weak"){ totalDofs = nprocs * ndofs; }

	// create and prepate 3D problem
	// 3D: tetrahedras and 1 quantity of interessed is assumed
	std::shared_ptr<dolfin::Expression> D3 = std::make_shared<TensorConstant>(rank, 0.013);
	std::shared_ptr<dolfin::Expression> init3D = std::make_shared<InitializerSphere>(0.5, 0.5, 0.5, 0.1, 1);
	dolfin::Timer t1("AAA create mesh 3D");
	int nxyz = pow((1.0 * totalDofs / 1.30) / 6, 1.0/3.0);
	std::shared_ptr<dolfin::Mesh> mesh3D = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(nxyz, nxyz, nxyz));
	t1.stop();
	std::shared_ptr<FisherProblem> problem3 = std::make_shared<FisherProblem>(rank, mesh3D, D3, 0.025, 0.0000001, 1);
	auto us3 = problem3->getUs();
	std::shared_ptr<dolfin::Function> u03 = us3.at(0);
	std::shared_ptr<dolfin::Function> u3 = us3.at(1);
	*u03 = *init3D;
	*u3 = *init3D;

	// Print simulation summary to console
	if (rank == 0)
	{
		std::stringstream ss;
		std::ofstream file("_INFO-Test-Summary");
		ss << "----------------------------------------------------------------" << std::endl <<
		 "Test problem summary" << std::endl <<
		"  Problem type:   "   << (problem_type == 1 ? "all" : (problem_type == 2 ? "cg as linsolver" : "unknown")) << " (3D)" << std::endl <<
		"  Scaling type:   "   << scaling_type << std::endl <<
		"  Num processes:  "  << nprocs << std::endl <<
		"  Mesh elements:	" <<  mesh3D->num_cells() << " (3D)" << std::endl <<
		"  Total dof:      " << us3.at(0)->function_space()->dim() << " (3D)" << std::endl <<
		"  Average dof per rank: " << us3.at(0)->function_space()->dim()/dolfin::MPI::size(mesh3D->mpi_comm()) << " (3D)" << std::endl <<
		"----------------------------------------------------------------" << std::endl;
		std::cout << ss.str() << std::endl;
		file << ss.str();
		file.close();
	}

	// different solver - preconditioner pairs
	///////////////////////////////////////////////////
	// NOTE:
	//		no parallel support for ilu/ icc preconditioner
	//		no support for petsc linear solver
	//		richardson and jacobi doesn't converge
	//
	dolfin::set_log_level(dolflog);

	std::stringstream iters;
	iters << "solver-iterations-" << nprocs << ".txt";
	std::ofstream iterfile(iters.str(), std::ios_base::trunc);
	iterfile.close();



	// tpye 1: run all, track time and iterations
	if(type == 1){
		runNewton(rank, nprocs, u3, problem3,  "gmres", "petsc_amg");
		runNewton(rank, nprocs, u3, problem3,  "gmres", "jacobi");
		runNewton(rank, nprocs, u3, problem3,  "gmres", "hypre_amg");
		runNewton(rank, nprocs, u3, problem3,  "gmres", "hypre_euclid");
		runNewton(rank, nprocs, u3, problem3,  "cg", "petsc_amg");
		runNewton(rank, nprocs, u3, problem3,  "cg", "jacobi");
		runNewton(rank, nprocs, u3, problem3,  "cg", "hypre_amg");
		runNewton(rank, nprocs, u3, problem3,  "cg", "hypre_euclid");
		runNewton(rank, nprocs, u3, problem3,  "richardson", "petsc_amg");
		runNewton(rank, nprocs, u3, problem3,  "richardson", "hypre_amg");
		runNewton(rank, nprocs, u3, problem3,  "richardson", "hypre_euclid");
		runNewton(rank, nprocs, u3, problem3,  "bicgstab", "hypre_amg");
		runNewton(rank, nprocs, u3, problem3,  "bicgstab", "hypre_euclid");
		runNewton(rank, nprocs, u3, problem3,  "minres", "hypre_amg");
		runNewton(rank, nprocs, u3, problem3,  "minres", "hypre_euclid");
		runNewton(rank, nprocs, u3, problem3,  "tfqmr", "hypre_amg");
		runNewton(rank, nprocs, u3, problem3,  "tfqmr", "hypre_euclid");
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
