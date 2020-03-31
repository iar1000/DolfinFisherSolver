#ifndef BRAIN_H
#define BRAIN_H

#include <vector>
#include <string>
#include <Eigen/Dense>

class Brain{

	int verbose_;
	std::string dataParentPath_;
	int dim_[3];
	std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> cms_;
	void readData();

public:
	// constructor for virtual brain
	// type:	0 = brainweb normal brain
	Brain(int verbose, std::string dataParentPath);
	// get concentration map and dimenstionality
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
		getConcentrationMap();
};

#endif
