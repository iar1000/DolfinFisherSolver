
#include <dolfin.h>
#include "../../DolfinFisherSolver/Tensors.h"
#include "../../DolfinFisherSolver/Initializers.h"
#include "../../DolfinFisherSolver/FisherProblem.h"
#include "../../DolfinFisherSolver/FisherNewtonContainer.h"


// automatic performance assestment of FisherSolver
int main(int argc, char* argv[]){
    // initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	dolfin::SubSystemsManager::init_mpi();
	dolfin::SubSystemsManager::init_petsc();

	// Default parameters
	dolfin::Parameters paras("application_parameters");
	paras.add("scaling_type", "weak", {"weak", "strong"});
	paras.add("ndofs", 500000);
	paras.add("output", true);

	// extract parameters
	const std::string problem_type = "fisher";
	const std::string scaling_type = paras["scaling_type"];
	const int ndofs = paras["ndofs"];
	const bool output = paras["output"];
	const std::string out_dir = "output/";

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
	std::shared_ptr<FisherProblem> problem2 = std::make_shared<FisherProblem>(rank, mesh2D, D2, 0.025, 0.0001, 1);
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
	std::shared_ptr<FisherProblem> problem3 = std::make_shared<FisherProblem>(rank, mesh3D, D3, 0.025, 0.0001, 1);
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

	std::shared_ptr<FisherNewtonContainer> fnc2 = std::make_shared<FisherNewtonContainer>(rank,
				mesh2D, init2D, D2, 0.025, 1, 0.00001);
	std::shared_ptr<FisherNewtonContainer> fnc3 = std::make_shared<FisherNewtonContainer>(rank,
					mesh3D, init3D, D3, 0.025, 1, 0.00001);

	// instanciate the different solvers
	std::vector<std::string> solver_methods = {"bicgstab", "cg", "gmres", "minres", "mumps", "petsc", "richardson", "superlu"};
	std::vector<std::string> solver_preconditioner = {"amg", "icc", "ilu", "jacobi", "petsc_amg", "sor"};
	for(int i = 0; i < solver_methods.size(); i++){
		for(int j = 0; j < solver_preconditioner.size(); j++){
			if(rank == 0){ std::cout << "solve " << solver_methods.at(i) << " + " << solver_preconditioner.at(j) << std::endl;}

			dolfin::Timer t2("AAA solve2D " + solver_methods.at(i) + "+" + solver_preconditioner.at(j));
			std::vector<std::pair<std::string, std::string>> paras;
			std::pair<std::string, std::string> sol = make_pair("linear_solver", solver_methods.at(i));
			std::pair<std::string, std::string> pre = make_pair("preconditioner", solver_preconditioner.at(j));
			paras.push_back(sol); paras.push_back(pre);
			fnc2->setSolverParameters(paras);
			fnc2->solve(0.1);
			t2.stop();

			dolfin::Timer t3("AAA solve3D " + solver_methods.at(i) + "+" + solver_preconditioner.at(j));
			paras.push_back(sol); paras.push_back(pre);
			fnc3->setSolverParameters(paras);
			fnc3->solve(0.1);
			t3.stop();
		}
	}



	//dolfin::list_linear_solver_methods();
	//dolfin::list_krylov_solver_preconditioners();


	 // Display timings
	 list_timings(dolfin::TimingClear::clear, {dolfin::TimingType::wall});

	 // output timings to text file
	 if(output){
		 std::ofstream timings;
		 timings.open("timings-latex-output", std::ios_base::trunc);
		 std::set<dolfin::TimingType> s = {dolfin::TimingType::wall};
		 timings << dolfin::timings(dolfin::TimingClear::keep, s).str_latex();
		 timings.close();
	 }

	MPI_Finalize();
}
