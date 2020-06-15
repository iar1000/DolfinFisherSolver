#ifndef BRAIN_H
#define BRAIN_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <dolfin.h>

#include "PrintableComponent.h"


class Brain : public PrintableComponent{

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

	// read in binary data from brainweb to local store
	void readBrainweb();
	// read in mesh data to local store
	void readMesh(std::shared_ptr<dolfin::Mesh> mesh);
	// compare virtual maps of mesh and brainweb datapoints
	std::vector<int> compareTranslated(std::shared_ptr<dolfin::Mesh> mesh, int translation[3], bool print, std::string print_name);
	// print out the given virtual map as slices
	void printVectorMatrix(std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> data, std::string name);
	// find points where white or grey tissue is located, store point list as csv
	void printBooleanSliceForAlphaShape();
public:
	// constructor for virtual brain
	Brain(int rank, int verbose, std::string dataParentPath);
	// find optimal translation within eps in every direction of the bounding translation
	// returns total tries, #misses, translations in x, y and z
	// print: create folder with hit-data of best translation for plotting with python script "plot_brain.py"
	std::vector<int> greedyOptimalTranslation(std::shared_ptr<dolfin::Mesh> mesh, int eps, bool print, std::string print_name);
	// evaluate fit of test meshes on the tissue concentration map
	// used to experiment on mesh
	void fitMesh();
	// get concentration map and dimenstionality
	// if slice between 0 and 180 returns matrix of this slice from brainweb
	// else: return all slices
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
		getConcentrationMap(int slice);
	// @override virtuals from PrintableComponent
	std::string asString();
};

#endif
