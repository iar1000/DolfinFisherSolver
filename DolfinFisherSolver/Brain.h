#ifndef BRAIN_H
#define BRAIN_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <dolfin.h>

class Brain{

	int rank_;							// rank of caller
	int verbose_;						// verbosity
	std::string dataParentPath_;		// path to brainweb data folder
	int dim_[3] = {181, 217, 181};		// dimensions of bounding box
	// brainweb
	std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> cms_;	// virtual brainweb concentration data maps
	double max_cm_[3];					// maximum values of dimensions x,y,z
	double min_cm_[3];					// minimum values of dimensions x,y,z
	double lengths_[3];					// length of span between extremas
	double middle_x;				// cutoff for left brainhalf
	// mesh
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> mms_;	// virtual mesh grid map
	double max_mesh_[3];					// maximum values of dimensions x,y,z
	double min_mesh_[3];					// minimum values of dimensions x,y,z
	double lengths_mesh_[3];				// length of span between extremas
	double middle_x_mesh_;				// cutoff for left brainhalf

	// read in binary data from brainweb
	void readBrainweb();
	// read in mesh data
	void readMesh(std::shared_ptr<dolfin::Mesh> mesh);
	// compare virtual maps
	std::vector<int> compareTranslated(std::shared_ptr<dolfin::Mesh> mesh, int translation[3], bool print);
	// print out the virtual maps
	void printVectorMatrix(std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> data, std::string name);

public:
	// constructor for virtual brain
	Brain(int rank, int verbose, std::string dataParentPath);
	// find optimal translation within eps of the bounding traslation
	// returns total tries, #misses, translations in x, y and z
	// print: create folder with hit-data of best translation for plotting with python script "plot_brain.py"
	std::vector<int> greedyOptimalTranslation(std::shared_ptr<dolfin::Mesh> mesh, int eps, bool print);
	// fit the tissue concentration map on the different meshes and fit best one
	// used to experiment on mesh
	void fitMesh();
	// get concentration map and dimenstionality
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
		getConcentrationMap();
};

#endif
