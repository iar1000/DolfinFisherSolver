#include <iostream>

#include "ValueMapper.h"

using namespace std;

ValueMapper::ValueMapper(int trank, int tdim_x, int tdim_y)
{
	rank = trank;
	verbose = 0;

	dim_x = tdim_x;
	dim_y = tdim_y;

	if(rank == 0){
		cout << "create ValueMapper: "
			"\n    verbose= 0" <<
			"\n    dimension X= " << to_string(dim_x) <<
			"\n    dimension Y= " << to_string(dim_y) << endl;
	};

};

// set verbosity on for debuging purposes
void ValueMapper::verbose_on()
{
	if (rank == 0) {cout << "ValueMapper::verbose_on: \n    set verbose= 1" << endl; };
	verbose = 1;
}

// create test value map for rect-100on100-res100.h test mesh
vector<Matrix<double, Dynamic, Dynamic>> ValueMapper::get_101on101_test_vm()
{
	vector<Matrix<double, Dynamic, Dynamic>> vm;

	// init valuemaps to form pattern
	Matrix<double, Dynamic, Dynamic> vm_w = MatrixXd::Ones(101, 101);
	vm_w.block<40, 101>(7, 0) = MatrixXd::Zero(40, 101);
	vm_w.block<40, 101>(54, 0) = MatrixXd::Zero(40, 101);

	Matrix<double, Dynamic, Dynamic> vm_g = MatrixXd::Zero(101, 101);
	vm_g.block<40, 101>(7, 0) = MatrixXd::Ones(40, 101);
	vm_g.block<40, 101>(54, 0) = MatrixXd::Ones(40, 101);

	vm.push_back(vm_w);
	vm.push_back(vm_g);
	return vm;
};
