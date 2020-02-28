
#include <iostream>
#include <Eigen/Dense>
#include <dolfin.h>

#include "../DolfinFisherSolver/Tensors.h"
#include "../DolfinFisherSolver/ReaderWriter.h"
#include "../DolfinFisherSolver/VariationalReactionDiffusion2D.h"
#include "../DolfinFisherSolver/Initializers.h"
#include "../DolfinFisherSolver/TimeStepper.h"
#include "../DolfinFisherSolver/ReactionDiffusionProblem.h"


// return a test concentraion map of size 100x100
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> test_cm_2D()
{
	// init valuemaps to form pattern
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vm_w = Eigen::MatrixXd::Ones(101, 101);
	vm_w.block<40, 101>(7, 0) = Eigen::MatrixXd::Zero(40, 101);
	vm_w.block<40, 101>(54, 0) = Eigen::MatrixXd::Zero(40, 101);

	return vm_w;
};

void testTensor(){
	std::cout << "Testing Tensor.cpp" << std::endl;
	// linear interpolation
	std::cout << "linear interpolation" << std::endl;
	double p1 = 0.3;
	double values1[2] = {0, 3};
	double should = 0.9;
	double is = linear_interpolation(p1, values1);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;
	double p11 = 0.75;
	double values11[2] = {3.8, 1.2};
	should = 1.85;
	is = linear_interpolation(p11, values11);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;

	// bilinear interpolation
	std::cout << "bilinear interpolation v1" << std::endl;
	double p2[2] = {0.3, 0.6};
	double values2[4] = {1, 3, 3, 1};
	should = 2.08;
	is = bilinear_interpolation(p2, values2);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;
	double p22[2] = {0.25, 0.7};
	double values22[4] = {0, 0, 4, 5};
	should = 1.175;
	is = bilinear_interpolation(p22, values22);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;
	double p222[2] = {0.43, 0.8};
	double values222[4] = {1, 7, 2, 12};
	should = 7.606;
	is = bilinear_interpolation(p222, values222);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;

	std::cout << "bilinear interpolation v2" << std::endl;
	should = 2.08;
	is = bilinear_interpolation_v2(p2, values2);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;
	should = 1.175;
	is = bilinear_interpolation_v2(p22, values22);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;

	// trilinear interpolation
	std::cout << "trilinear interpolation v1" << std::endl;
	double p3[3] = {0.2, 0.47, 0.9};
	double values3[8] = {1, 1, 1, 1, 3, 3, 3, 3};
	should = 2.8;
	is = trilinear_interpolation(p3, values3);
	std::cout << "	test result: " << (is == should ? "passed": "failed") << " (" << should << "," << is << ")" << std::endl;

	double res1 = trilinear_interpolation(p3, values3);
	double res2 = trilinear_interpolation_v2(p3, values3);
	std::cout << "compare trilinear interpolation v1 vs. v2" << std::endl;
	std::cout << "	test result: " << (res1 == res2 ? "equal result" : "different result") <<
			"  (" << res1 << "," << res2 << ")" << std::endl;
	double p33[3] = {0.72, 0.5, 0.12};
	double values33[8] = {1, 3, 3, 8, 1, 5, 0, 3};
	res1 = trilinear_interpolation(p33, values33);
	res2 = trilinear_interpolation_v2(p33, values33);

	std::cout << "	test result: " << (res1 == res2 ? "equal result" : "different result") <<
			"  (" << res1 << "," << res2 << ")" << std::endl;

}

void testReaderWriter(){
	std::cout << "Testing ReaderWriter.cpp" << std::endl;

	// string split
	std::cout << "splitString" << std::endl;
	std::string toSplit = "should.should-not.should_not be-split.should";
	std::vector<std::string> tokens = splitString(toSplit, '.');
	std::cout << "	result of " << toSplit << std::endl;
	for(unsigned i = 0; i < tokens.size(); i++){
		std::cout << "	  " << tokens.at(i) << std::endl;
	}
	std::cout << " " << std::endl;

	ReaderWriter rw = ReaderWriter(0, "output", "../mesh");

	// load mesh
	std::cout << "loadMesh" << std::endl;
	std::shared_ptr<dolfin::Mesh> meshWrong = std::make_shared<dolfin::Mesh>();
	auto wrong = rw.loadMesh(meshWrong, "wrong.h5");
	std::cout << "	test loading wrong.h5: " << (wrong.first ? "failed" : "passed") << std::endl;
	wrong = rw.loadMesh(meshWrong, "format.h5.xyz");
	std::cout << "	test loading format.h5.xyz: " << (wrong.first ? "failed" : "passed") << std::endl;
	wrong = rw.loadMesh(meshWrong, "format.xyz");
	std::cout << "	test loading format.xyz: " << (wrong.first ? "failed" : "passed") << std::endl;

	std::shared_ptr<dolfin::Mesh> mesh1 = std::make_shared<dolfin::Mesh>();
	auto right = rw.loadMesh(mesh1, "rect-100on100-res-100.h5");
	std::cout << "	test loading rect-100on100-res-100.h5: " << (right.first ? "passed" : "failed") <<
			", dimensions (" << right.second << "," << mesh1->geometry().dim() << ")" << std::endl;
	std::shared_ptr<dolfin::Mesh> mesh2 = std::make_shared<dolfin::Mesh>();
	right = rw.loadMesh(mesh2, "box-10on10on10-res-50.h5");
	std::cout << "	test loading box-10on10on10-res-50.h5: " << (right.first ? "passed" : "failed") <<
			", dimension (" << right.second << "," << mesh2->geometry().dim() << ")" << std::endl;
	std::cout << std::endl;

}

void testMeshFunctions(){
	ReaderWriter putput = ReaderWriter(0, "output", "../mesh");

	// mesh read-in
	std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>();
	std::pair<bool, int> meshInfo = putput.loadMesh(mesh, "rect-100on100-res-100.h5");
	auto cellType = mesh->type().cell_type();
	std::cout << "mesh: " << std::endl <<
			"	cell type = "  << dolfin::CellType::type2string(cellType) << std::endl <<
			"	cell description = "  << mesh->type().description(true) << std::endl <<
			"	cell top dim = "  << mesh->type().dim() << std::endl <<
			"	cell num verticies = "  << mesh->type().num_vertices() << std::endl <<
			"	mesh num cells (local) = "  << mesh->num_cells() << std::endl <<
			"	mesh num edges (local) = "  << mesh->num_edges() << std::endl <<
			"	mesh num faces (local) = "  << mesh->num_faces() << std::endl;


	std::cout << std::endl;
	auto V = std::make_shared<VariationalReactionDiffusion2D::FunctionSpace>(mesh);
	auto dofs = V->dofmap();
	std::cout << "dofmap: " << std::endl <<
			"	global dimension = "  << dofs->global_dimension() << std::endl <<
			"	cell 10 dimension = "  << dofs->cell_dimension(10) << std::endl <<
			"	max cell dimension = "  << dofs->max_cell_dimension() << std::endl <<
			"	num facet dofs = "  << dofs->num_facet_dofs() << std::endl <<
			"	size dofs list = "  << dofs->dofs().size() << std::endl;


}

void testAdaptiveMesh(){

	// initialize MPI
	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>();
	auto hdf5 = dolfin::HDF5File(MPI_COMM_WORLD, "../mesh/rect-10on10-res-25.h5", std::string("r"));
	hdf5.read(*mesh, "/mesh", false);

	// create initial condition
	std::shared_ptr<InitializerCircle> initialCondition;
	initialCondition = std::make_shared<InitializerCircle>(5, 5, 1, 1);
	if(rank == 0){ std::cout << initialCondition->asString() << std::endl; };

	std::shared_ptr<TensorConstant> DConstant = std::make_shared<TensorConstant>(rank, 0.5);
	if(rank == 0){ std::cout << DConstant->asString() << std::endl; };

	std::shared_ptr<ReactionDiffusionProblem> problem = std::make_shared<ReactionDiffusionProblem>(rank, mesh, DConstant, 0.001, 0.001, 1);
	if(rank == 0){ std::cout << problem->asString() << std::endl;};

	std::shared_ptr<dolfin::NewtonSolver> solver = std::make_shared<dolfin::NewtonSolver>();
	solver->parameters["linear_solver"] = "lu";
	solver->parameters["convergence_criterion"] = "incremental";
	solver->parameters["maximum_iterations"] = 50;
	solver->parameters["relative_tolerance"] = 1e-10;
	solver->parameters["absolute_tolerance"] = 1e-10;

	if(rank == 0){
		std::cout << solver->parameters.str(true) << std::endl;
	}

	return;

	TimeStepper timeStepper = TimeStepper(rank, problem, solver, 10, 0.00000001, 0.01, 0.001, 1);
	if(rank == 0){ std::cout << timeStepper.asString() << std::endl; };

	// create output files
	std::shared_ptr<dolfin::File> file = std::make_shared<dolfin::File>("output/adaptiveMesh/out.pvd");

	// run simulation
	RuntimeTracker tracker3 = timeStepper.run(1, 3, initialCondition,
			file, "output/adaptiveMesh/iterationdata.csv",
			-1, 0.001);

	//AdaptiveLinearVariationalSolver solver(problem, M);


	MPI_Finalize(); //seems to trigger an abort
}

int main()
{
	testAdaptiveMesh();

};
