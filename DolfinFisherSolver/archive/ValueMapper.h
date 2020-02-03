#ifndef VALUEMAPPER_H
#define VALUEMAPPER_H

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

class ValueMapper
{
	int verbose, rank;
	int dim_x, dim_y;
public:
	ValueMapper(int rank, int tdim_x, int tdim_y);
	void verbose_on(); // set verbosity on for debugging purposes

	std::vector<Matrix<double, Dynamic, Dynamic>> get_101on101_test_vm();
};

#endif
