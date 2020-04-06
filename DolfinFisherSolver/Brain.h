#ifndef BRAIN_H
#define BRAIN_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <dolfin.h>

class Brain{

	int rank_;							// rank of caller
	int verbose_;						// verbosity
	// brainweb raw specific
	std::string dataParentPath_;		// path to brainweb data folder
	int dim_[3] = {181, 217, 181};		// dimensions of brainweb data
	std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> cms_;	// virtual brainweb concentration data maps
	// brainweb additional info
	double max_cm_[3];					// maximum values of dimensions x,y,z
	double min_cm_[3];					// minimum values of dimensions x,y,z
	double lengths_[3];					// length of span between extremas
	// read in binary data from brainweb
	void readData();

public:
	// constructor for virtual brain
	Brain(int rank, int verbose, std::string dataParentPath);
	// fit the tissue concentration map on the mesh version 1 for 3D
	void fitMesh(std::shared_ptr<dolfin::Mesh> mesh);
	// get concentration map and dimenstionality
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
		getConcentrationMap();
};

#endif
