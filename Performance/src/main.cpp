
#include <dolfin.h>
#include <fstream>
#include "../../DolfinFisherSolver/Tensors.h"
#include "../../DolfinFisherSolver/Initializers.h"
#include "../../DolfinFisherSolver/FisherProblem.h"
#include "../../DolfinFisherSolver/FisherNewtonContainer.h"


// automatic performance assestment of FisherSolver
int main(int argc, char* argv[]){
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	if(rank==0){
		dolfin::list_krylov_solver_preconditioners();
		dolfin::list_linear_solver_methods();
	}

	// Parse command line options (will intialise PETSc if any PETSc
	// options are present, e.g. --petsc.pc_type=jacobi)
	dolfin::parameters.parse(argc, argv);
	dolfin::SubSystemsManager::init_petsc();

	// Default parameters
	dolfin::Parameters application_parameters("application_parameters");
	application_parameters.add("scaling_type", "weak", {"weak", "strong"});
	application_parameters.add("ndofs", 500000);
	application_parameters.add("output", true);

	// Update from command line
	application_parameters.parse(argc, argv);

	// extract parameters
	const std::string problem_type = "fisher";
	const std::string scaling_type = application_parameters["scaling_type"];
	const int ndofs = application_parameters["ndofs"];
	const bool output = application_parameters["output"];
	const std::string out_dir = "output/";

	 // Set mesh partitioner
	 dolfin::parameters["mesh_partitioner"] = "SCOTCH";

	// create constant diffusion tensors
	std::shared_ptr<dolfin::Expression> D2 = std::make_shared<TensorConstant>(rank, 0.013);
	std::shared_ptr<dolfin::Expression> D3 = std::make_shared<TensorConstant>(rank, 0.013);

	// create inital conditions
	std::shared_ptr<dolfin::Expression> init2D = std::make_shared<InitializerCircle>(0.5, 0.5, 0.1, 1);
	std::shared_ptr<dolfin::Expression> init3D = std::make_shared<InitializerSphere>(0.5, 0.5, 0.5, 0.1, 1);

	// calculate dofs
	int totalDofs = ndofs;
	if(scaling_type == "weak"){ totalDofs = nprocs * ndofs; }

	// create and prepare 2D problem
	// 2D: triangles and 1 quantity of interessed is assumed
	dolfin::Timer t0("AAA create mesh 2D");
	int nxy = sqrt((totalDofs / 2) / 2);
	std::shared_ptr<dolfin::Mesh> mesh2D = std::make_shared<dolfin::Mesh>(dolfin::UnitSquareMesh(nxy, nxy, "right"));
	std::shared_ptr<FisherProblem> problem2 = std::make_shared<FisherProblem>(rank, mesh2D, D2, 0.025, 0.0000001, 1);
	auto us2 = problem2->getUs();
	std::shared_ptr<dolfin::Function> u02 = us2.at(0);
	std::shared_ptr<dolfin::Function> u2 = us2.at(1);
	*u02 = *init2D;
	*u2 = *init2D;
	t0.stop();

	// create and prepate 3D problem
	// 3D: tetrahedras and 1 quantity of interessed is assumed
	dolfin::Timer t1("AAA create mesh 3D");
	int nxyz = pow((1.0 * totalDofs / 1.30) / 6, 1.0/3.0);
	std::shared_ptr<dolfin::Mesh> mesh3D = std::make_shared<dolfin::Mesh>(dolfin::UnitCubeMesh(nxyz, nxyz, nxyz));
	std::shared_ptr<FisherProblem> problem3 = std::make_shared<FisherProblem>(rank, mesh3D, D3, 0.025, 0.0000001, 1);
	auto us3 = problem3->getUs();
	std::shared_ptr<dolfin::Function> u03 = us3.at(0);
	std::shared_ptr<dolfin::Function> u3 = us3.at(1);
	*u03 = *init3D;
	*u3 = *init3D;
	t1.stop();

	// Print simulation summary to console
	if (rank == 0)
	{
		std::cout << "----------------------------------------------------------------" << std::endl;
		std::cout << "Test problem summary" << std::endl;
		std::cout << "  Problem type:   "   << problem_type << std::endl;
		std::cout << "  Scaling type:   "   << scaling_type << std::endl;
		std::cout << "  Num processes:  "  << nprocs << std::endl;
		std::cout << "  Mesh elements:	" << mesh2D->num_cells() << " (2D)/ " << mesh3D->num_cells() << " (3D)" << std::endl;
		std::cout << "  Total dof:      " <<  us2.at(0)->function_space()->dim()<< " (2D)/ " << us3.at(0)->function_space()->dim() << " (3D)" << std::endl;
		std::cout << "  Average dof per rank: " << us2.at(0)->function_space()->dim()/dolfin::MPI::size(mesh2D->mpi_comm()) << " (2D)/ " <<
				us3.at(0)->function_space()->dim()/dolfin::MPI::size(mesh3D->mpi_comm()) << " (3D)" << std::endl;
		std::cout << "----------------------------------------------------------------" << std::endl;
	}

	// different solver - preconditioner pairs
	///////////////////////////////////////////////////
	dolfin::set_log_level(30);
	// gmres
	if(rank==0){ std::cout << "BBB solve gmres + pets_amg" << std::endl;}
	dolfin::Timer t2("BBB solve gmres + pets_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver2 = std::make_shared<dolfin::NewtonSolver>();
	solver2->parameters["convergence_criterion"] = "incremental";
	solver2->parameters["linear_solver"] = "gmres";
	solver2->parameters["preconditioner"] = "petsc_amg";
	solver2->solve(*problem2, *u2->vector());
	t2.stop();

	if(rank==0){ std::cout << "BBB solve gmres + jacobi" << std::endl;}
	dolfin::Timer t3("BBB solve gmres + jacobi");
	std::shared_ptr<dolfin::NewtonSolver> solver3 = std::make_shared<dolfin::NewtonSolver>();
	solver3->parameters["convergence_criterion"] = "incremental";
	solver3->parameters["linear_solver"] = "gmres";
	solver3->parameters["preconditioner"] = "jacobi";
	solver3->solve(*problem2, *u2->vector());
	t3.stop();

	/* No parallel ilu preconditioner support
	if(rank==0){ std::cout << "BBB solve gmres + ilu" << std::endl;}
	dolfin::Timer t4("BBB solve gmres + ilu");
	std::shared_ptr<dolfin::NewtonSolver> solver4 = std::make_shared<dolfin::NewtonSolver>();
	solver4->parameters["convergence_criterion"] = "incremental";
	solver4->parameters["linear_solver"] = "gmres";
	solver4->parameters["preconditioner"] = "ilu";
	solver4->solve(*problem2, *u2->vector());
	t4.stop();

	if(rank==0){ std::cout << "BBB solve gmres + icc" << std::endl;}
	dolfin::Timer t5("BBB solve gmres + icc");
	std::shared_ptr<dolfin::NewtonSolver> solver5 = std::make_shared<dolfin::NewtonSolver>();
	solver5->parameters["convergence_criterion"] = "incremental";
	solver5->parameters["linear_solver"] = "gmres";
	solver5->parameters["preconditioner"] = "icc";
	solver5->solve(*problem2, *u2->vector());
	t5.stop();
	 */

	if(rank==0){ std::cout << "BBB solve gmres + hypre_amg" << std::endl;}
	dolfin::Timer t6("BBB solve gmres + hypre_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver6 = std::make_shared<dolfin::NewtonSolver>();
	solver6->parameters["convergence_criterion"] = "incremental";
	solver6->parameters["linear_solver"] = "gmres";
	solver6->parameters["preconditioner"] = "hypre_amg";
	solver6->solve(*problem2, *u2->vector());
	t6.stop();

	if(rank==0){ std::cout << "BBB solve gmres + hypre_euclid" << std::endl;}
	dolfin::Timer t7("BBB solve gmres + hypre_euclid");
	std::shared_ptr<dolfin::NewtonSolver> solver7 = std::make_shared<dolfin::NewtonSolver>();
	solver7->parameters["convergence_criterion"] = "incremental";
	solver7->parameters["linear_solver"] = "gmres";
	solver7->parameters["preconditioner"] = "hypre_euclid";
	solver7->solve(*problem2, *u2->vector());
	t7.stop();

	// cg
	if(rank==0){ std::cout << "BBB solve cg + pets_amg" << std::endl;}
	dolfin::Timer t8("BBB solve cg + pets_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver8 = std::make_shared<dolfin::NewtonSolver>();
	solver8->parameters["convergence_criterion"] = "incremental";
	solver8->parameters["linear_solver"] = "cg";
	solver8->parameters["preconditioner"] = "petsc_amg";
	solver8->solve(*problem2, *u2->vector());
	t8.stop();

	if(rank==0){ std::cout << "BBB solve cg + jacobi" << std::endl;}
	dolfin::Timer t9("BBB solve cg + jacobi");
	std::shared_ptr<dolfin::NewtonSolver> solver9 = std::make_shared<dolfin::NewtonSolver>();
	solver9->parameters["convergence_criterion"] = "incremental";
	solver9->parameters["linear_solver"] = "cg";
	solver9->parameters["preconditioner"] = "jacobi";
	solver9->solve(*problem2, *u2->vector());
	t9.stop();

	/* No parallel ilu preconditioner support
	if(rank==0){ std::cout << "BBB solve cg + ilu" << std::endl;}
	dolfin::Timer t10("BBB solve cg + ilu");
	std::shared_ptr<dolfin::NewtonSolver> solver10 = std::make_shared<dolfin::NewtonSolver>();
	solver10->parameters["convergence_criterion"] = "incremental";
	solver10->parameters["linear_solver"] = "cg";
	solver10->parameters["preconditioner"] = "ilu";
	solver10->solve(*problem2, *u2->vector());
	t10.stop();

	if(rank==0){ std::cout << "BBB solve cg + icc" << std::endl;}
	dolfin::Timer t11("BBB solve cg + icc");
	std::shared_ptr<dolfin::NewtonSolver> solver11 = std::make_shared<dolfin::NewtonSolver>();
	solver11->parameters["convergence_criterion"] = "incremental";
	solver11->parameters["linear_solver"] = "cg";
	solver11->parameters["preconditioner"] = "icc";
	solver11->solve(*problem2, *u2->vector());
	t11.stop();
	*/

	if(rank==0){ std::cout << "BBB solve cg + hypre_amgg" << std::endl;}
	dolfin::Timer t12("BBB solve cg + hypre_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver12 = std::make_shared<dolfin::NewtonSolver>();
	solver12->parameters["convergence_criterion"] = "incremental";
	solver12->parameters["linear_solver"] = "cg";
	solver12->parameters["preconditioner"] = "hypre_amg";
	solver12->solve(*problem2, *u2->vector());
	t12.stop();

	if(rank==0){ std::cout << "BBB solve cg + hypre_euclid" << std::endl;}
	dolfin::Timer t13("BBB solve cg + hypre_euclid");
	std::shared_ptr<dolfin::NewtonSolver> solver13 = std::make_shared<dolfin::NewtonSolver>();
	solver13->parameters["convergence_criterion"] = "incremental";
	solver13->parameters["linear_solver"] = "cg";
	solver13->parameters["preconditioner"] = "hypre_euclid";
	solver13->solve(*problem2, *u2->vector());
	t13.stop();

	// petsc

	if(rank==0){ std::cout << "BBB solve petsc" << std::endl;}
	dolfin::Timer t14("BBB solve petsc");
	std::shared_ptr<dolfin::NewtonSolver> solver14 = std::make_shared<dolfin::NewtonSolver>();
	solver14->parameters["convergence_criterion"] = "incremental";
	solver14->parameters["linear_solver"] = "petsc";
	solver14->solve(*problem2, *u2->vector());
	t14.stop();

	/* not possible to use preconditioner with this solver
	if(rank==0){ std::cout << "BBB solve petsc + jacobi" << std::endl;}
	dolfin::Timer t15("BBB solve petsc + jacobi");
	std::shared_ptr<dolfin::NewtonSolver> solver15 = std::make_shared<dolfin::NewtonSolver>();
	solver15->parameters["convergence_criterion"] = "incremental";
	solver15->parameters["linear_solver"] = "petsc";
	solver15->parameters["preconditioner"] = "jacobi";
	solver15->solve(*problem2, *u2->vector());
	t15.stop();


	if(rank==0){ std::cout << "BBB solve petsc + ilu" << std::endl;}
	dolfin::Timer t16("BBB solve petsc + ilu");
	std::shared_ptr<dolfin::NewtonSolver> solver16 = std::make_shared<dolfin::NewtonSolver>();
	solver16->parameters["convergence_criterion"] = "incremental";
	solver16->parameters["linear_solver"] = "petsc";
	solver16->parameters["preconditioner"] = "ilu";
	solver16->solve(*problem2, *u2->vector());
	t16.stop();

	if(rank==0){ std::cout << "BBB solve petsc + icc" << std::endl;}
	dolfin::Timer t17("BBB solve petsc + icc");
	std::shared_ptr<dolfin::NewtonSolver> solver17 = std::make_shared<dolfin::NewtonSolver>();
	solver17->parameters["convergence_criterion"] = "incremental";
	solver17->parameters["linear_solver"] = "petsc";
	solver17->parameters["preconditioner"] = "icc";
	solver17->solve(*problem2, *u2->vector());
	t17.stop();

	if(rank==0){ std::cout << "BBB solve petsc + hypre_amg" << std::endl;}
	dolfin::Timer t18("BBB solve petsc + hypre_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver18 = std::make_shared<dolfin::NewtonSolver>();
	solver18->parameters["convergence_criterion"] = "incremental";
	solver18->parameters["linear_solver"] = "petsc";
	solver18->parameters["preconditioner"] = "hypre_amg";
	solver18->solve(*problem2, *u2->vector());
	t18.stop();

	if(rank==0){ std::cout << "BBB solve petsc + hypre_euclid" << std::endl;}
	dolfin::Timer t19("BBB solve petsc + hypre_euclid");
	std::shared_ptr<dolfin::NewtonSolver> solver19 = std::make_shared<dolfin::NewtonSolver>();
	solver19->parameters["convergence_criterion"] = "incremental";
	solver19->parameters["linear_solver"] = "petsc";
	solver19->parameters["preconditioner"] = "hypre_euclid";
	solver19->solve(*problem2, *u2->vector());
	t19.stop();
	 */

	// richardson
	if(rank==0){ std::cout << "BBB solve richardson + pets_amg" << std::endl;}
	dolfin::Timer t20("BBB solve richardson + pets_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver20 = std::make_shared<dolfin::NewtonSolver>();
	solver20->parameters["convergence_criterion"] = "incremental";
	solver20->parameters["linear_solver"] = "richardson";
	solver20->parameters["preconditioner"] = "petsc_amg";
	solver20->solve(*problem2, *u2->vector());
	t20.stop();

	/* doesn't converge
	if(rank==0){ std::cout << "BBB solve richardson + jacobi" << std::endl;}
	dolfin::Timer t21("BBB solve richardson + jacobi");
	std::shared_ptr<dolfin::NewtonSolver> solver21 = std::make_shared<dolfin::NewtonSolver>();
	solver21->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	solver21->parameters["convergence_criterion"] = "incremental";
	solver21->parameters["linear_solver"] = "richardson";
	solver21->parameters["preconditioner"] = "jacobi";
	solver21->solve(*problem2, *u2->vector());
	t21.stop();
	*/

	/* No parallel ilu preconditioner support
	if(rank==0){ std::cout << "BBB solve richardson + ilu" << std::endl;}
	dolfin::Timer t22("BBB solve richardson + ilu");
	std::shared_ptr<dolfin::NewtonSolver> solver22 = std::make_shared<dolfin::NewtonSolver>();
	solver22->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	solver22->parameters["convergence_criterion"] = "incremental";
	solver22->parameters["linear_solver"] = "richardson";
	solver22->parameters["preconditioner"] = "ilu";
	solver22->solve(*problem2, *u2->vector());
	t22.stop();

	if(rank==0){ std::cout << "BBB solve richardson + icc" << std::endl;}
	dolfin::Timer t23("BBB solve richardson + icc");
	std::shared_ptr<dolfin::NewtonSolver> solver23 = std::make_shared<dolfin::NewtonSolver>();
	solver23->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	solver23->parameters["convergence_criterion"] = "incremental";
	solver23->parameters["linear_solver"] = "richardson";
	solver23->parameters["preconditioner"] = "icc";
	solver23->solve(*problem2, *u2->vector());
	t23.stop();
	*/

	if(rank==0){ std::cout << "BBB solve richardson + hypre_amg" << std::endl;}
	dolfin::Timer t24("BBB solve richardson + hypre_amg");
	std::shared_ptr<dolfin::NewtonSolver> solver24 = std::make_shared<dolfin::NewtonSolver>();
	solver24->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	solver24->parameters["convergence_criterion"] = "incremental";
	solver24->parameters["linear_solver"] = "richardson";
	solver24->parameters["preconditioner"] = "hypre_amg";
	solver24->solve(*problem2, *u2->vector());
	t24.stop();

	if(rank==0){ std::cout << "BBB solve richardson + hypre_euclid" << std::endl;}
	dolfin::Timer t25("BBB solve richardson + hypre_euclid");
	std::shared_ptr<dolfin::NewtonSolver> solver25 = std::make_shared<dolfin::NewtonSolver>();
	solver25->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	solver25->parameters["convergence_criterion"] = "incremental";
	solver25->parameters["linear_solver"] = "richardson";
	solver25->parameters["preconditioner"] = "hypre_euclid";
	solver25->solve(*problem2, *u2->vector());
	t25.stop();



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

		 dolfin::dump_timings_to_xml("timings.xml", dolfin::TimingClear::keep);
	 }

	MPI_Finalize();
}
