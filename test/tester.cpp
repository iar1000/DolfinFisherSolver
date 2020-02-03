
#include <iostream>
#include <Eigen/Dense>
#include <dolfin.h>

#include "../DolfinFisherSolver/Tensors.h"
#include "../DolfinFisherSolver/ReaderWriter.h"

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

	// load file
	std::cout << "loadFile" << std::endl;
	std::shared_ptr<dolfin::File> file1;
	bool fileLoaded = rw.loadFile(file1, "should-work", "out", "pvd");
	std::cout << "	test loading should-work/out.pvd " << (fileLoaded ? "success" : "failed") << std::endl;
	std::shared_ptr<dolfin::File> file11;
	fileLoaded = rw.loadFile(file11, "should-work", "out2", "pvd");
	std::cout << "	test loading should-work/out2.pvd " << (fileLoaded ? "success" : "failed") << std::endl;

	std::shared_ptr<dolfin::File> file2;
	fileLoaded = rw.loadFile(file2, "", "out", "pvd");
	std::cout << "	test loading /out.pvd " << (fileLoaded ? "failed" : "success") << std::endl;
	std::shared_ptr<dolfin::File> file3;
	fileLoaded = rw.loadFile(file3, "should-not-work", "out", "bla");
	std::cout << "	test loading should-not-work/out.bla " << (fileLoaded ? "failed" : "success") << std::endl;

}

int main()
{
	testReaderWriter();


};
