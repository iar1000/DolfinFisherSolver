
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <sys/stat.h>
#include <sstream>


#include "Brain.h"

Brain::Brain(int rank, int verbose, std::string dataParentPath){
	rank_ = rank;
	verbose_ = verbose;
	dataParentPath_ = dataParentPath;
	// initialize extremas and read in data
	middle_x = 88;
	for(int i = 0; i < 3; i++){
		max_cm_[i] = -99999;
		min_cm_[i] = 999999;
		max_mesh_[i] = -99999;
		min_mesh_[i] = 999999;
	}
	readBrainweb();
}

std::string Brain::asString(){
	std::stringstream ss3;
	ss3 << "Virtual Brain: " << std::endl << "\t Concentration Map size= " << dim_[0] << "x" << dim_[1] << "x" << dim_[2] << std::endl <<
				"\t Read-in border points of CM= " << max_cm_[0] << ", " << max_cm_[1] << ", " <<  max_cm_[2] << " (max x,y,z), " <<
				 min_cm_[0] << ", " << min_cm_[1] << ", " <<  min_cm_[2] << " (min x,y,z)" << std::endl;
	return ss3.str();
}



void Brain::readBrainweb(){
	if(rank_ == 0 && verbose_ > 3){ std::cout << "	read in BrainWeb tissue data..." << std::endl; }
	std::string greyPath = dataParentPath_ + "/phantom_1.0mm_normal_gry.rawb";
	std::string whitePath = dataParentPath_ + "/phantom_1.0mm_normal_wht.rawb";
	std::string paths[2] = {whitePath, greyPath};
	double range = 255;
	int x = dim_[0];
	int y = dim_[1];
	int z = dim_[2];
	int plane = x*y;
	char *buf = new char[1];

	for(int s = 0; s < 2; s++){
		// create empty concentration map
		std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> worker;
		for(int k = 0; k < z; k++){
			Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> slice = Eigen::MatrixXd::Zero(y, x);
			worker.push_back(slice);
		}
		// read in tissue data
		std::ifstream data(paths[s], std::ios::binary);
		data.read(buf, 1);
		int i = 0;
		int slice_counter = 0;
		int z_curr = 0;
		while(data){
			// information about current voxel
			double x_curr = slice_counter%x;
			double y_curr = slice_counter/x;
			double value = ((double)((unsigned char)*buf)) / range;
			// update value map
			worker.at(z_curr)(y_curr, x_curr) = value;
			i++;
			slice_counter++;
			if(slice_counter == plane){
				slice_counter = 0;
				z_curr++;
			}
			data.read(buf, 1);
			// update additional information about extremas
			if(value > 0){
				if(x_curr > max_cm_[0]){ max_cm_[0] = x_curr; }
				if(x_curr < min_cm_[0]){ min_cm_[0] = x_curr; }
				if(y_curr > max_cm_[1]){ max_cm_[1] = y_curr; }
				if(y_curr < min_cm_[1]){ min_cm_[1] = y_curr; }
				if(z_curr > max_cm_[2]){ max_cm_[2] = z_curr; }
				if(z_curr < min_cm_[2]){ min_cm_[2] = z_curr; }
			}
		}
		data.close();
		cms_.push_back(worker);
		if(rank_ == 0 && verbose_ > 3){ std::cout << "		" << i << " data points from " << paths[s] << std::endl; }
	}
	// update span lengths
	for(int i = 0; i < 3; i++){ lengths_[i] = max_cm_[i] - min_cm_[i]; };
	if(rank_ == 0 && verbose_ > 3){
		std::cout << "	tissue bounds:" << std::endl <<
				"		x: [" << min_cm_[0] << ", " << max_cm_[0] << "]" << " (" << lengths_[0] << ")" << std::endl <<
				"		y: [" << min_cm_[1] << ", " << max_cm_[1] << "]" << " (" << lengths_[1] << ")" << std::endl <<
				"		z: [" << min_cm_[2] << ", " << max_cm_[2] << "]" << " (" << lengths_[2] << ")" << std::endl;
	}
}

void Brain::readMesh(std::shared_ptr<dolfin::Mesh> mesh){
	// setup of variables
	if(rank_ == 0 && verbose_ > 3){ std::cout << "	read in mesh data..." << std::endl; }
	int x = dim_[0];
	int y = dim_[1];
	int z = dim_[2];
	char *buf = new char[1];

	// initialize local empty concentration map
	mms_.clear();
	for(int i = 0; i < 3; i++){
		max_mesh_[i] = -99999;
		min_mesh_[i] = 999999;
	}
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> worker;
	for(int k = 0; k < z; k++){
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> slice = Eigen::MatrixXd::Zero(y, x);
		worker.push_back(slice);
	}

	// go through coordinates
	auto coordinates = mesh->coordinates();
	for(int i = 0; i < coordinates.size(); i += 3){
		int x[2] = {(int)floor(coordinates[i]), (int)ceil(coordinates[i])};
		int y[2] = {(int)floor(coordinates[i+1]), (int)ceil(coordinates[i+1])};
		int z[2] = {(int)floor(coordinates[i+2]), (int)ceil(coordinates[i+2])};

		// set markers
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){	// virtual mesh grid map
				worker.at(z[j])(y[k], x[0]) = 1;
				worker.at(z[j])(y[k], x[1]) = 1;
			}
		}

		// update bounding box
		if(coordinates[i] > max_mesh_[0]){ max_mesh_[0] = coordinates[i]; }
		if(coordinates[i] < min_mesh_[0]){ min_mesh_[0] = coordinates[i]; }
		if(coordinates[i+1] > max_mesh_[1]){ max_mesh_[1] = coordinates[i+1]; }
		if(coordinates[i+1] < min_mesh_[1]){ min_mesh_[1] = coordinates[i+1]; }
		if(coordinates[i+2] > max_mesh_[2]){ max_mesh_[2] = coordinates[i+2]; }
		if(coordinates[i+2] < min_mesh_[2]){ min_mesh_[2] = coordinates[i+2]; }
	}
	mms_ = worker;

	// extract important parameters
	for(int i = 0; i < 3; i++){ lengths_mesh_[i] = max_mesh_[i] - min_mesh_[i]; }
	// print dimensions
	if(rank_ == 0 && verbose_ > 3){
		std::cout << "	mesh bounds: " << std::endl <<
				"		x: [" << min_mesh_[0] << ", " << max_mesh_[0] << "]" << " (" << lengths_mesh_[0] << ")" << std::endl <<
				"		y: [" << min_mesh_[1] << ", " << max_mesh_[1] << "]" << " (" << lengths_mesh_[1] << ")" << std::endl <<
				"		z: [" << min_mesh_[2] << ", " << max_mesh_[2] << "]" << " (" << lengths_mesh_[2] << ")" << std::endl;
	}
}

std::vector<int> Brain::compareTranslated(std::shared_ptr<dolfin::Mesh> mesh, int translation[3], bool print, std::string print_name){
	// create empty comparison map
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> worker;
	for(int k = 0; k < dim_[2]; k++){
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> slice = Eigen::MatrixXd::Zero(dim_[1], dim_[0]);
		worker.push_back(slice);
	}
	// counter variables
	int hit_white = 0;
	int hit_grey = 0;
	int miss = 0;
	int tries = 0;
	// coordinates
	auto coordinates = mesh->coordinates();
	for(int i = 0; i < coordinates.size(); i += 3){
		// get coordinates
		int x[2] = {(int)floor(coordinates[i]), (int)ceil(coordinates[i])};
		int y[2] = {(int)floor(coordinates[i+1]), (int)ceil(coordinates[i+1])};
		int z[2] = {(int)floor(coordinates[i+2]), (int)ceil(coordinates[i+2])};
		// translate
		for(int t = 0; t < 2; t++){
			x[t] += translation[0];
			y[t] += translation[1];
			z[t] += translation[2];
		}
		// set markers
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){
				tries += 2;
				// if grey or white hits then mark the spot
				if(cms_[0].at(z[j])(y[k], x[0]) != 0 ||
						cms_[1].at(z[j])(y[k], x[0]) != 0){
					worker.at(z[j])(y[k], x[0]) = 1;
					// counter
					if(cms_[0].at(z[j])(y[k], x[0]) != 0){ hit_white++; }
					if(cms_[1].at(z[j])(y[k], x[0]) != 0){ hit_grey++; }
				}
				else{ miss++; }
				if(cms_[0].at(z[j])(y[k], x[1]) != 0 ||
						cms_[1].at(z[j])(y[k], x[1]) != 0){
					worker.at(z[j])(y[k], x[1]) = 1;
					// counter
					if(cms_[0].at(z[j])(y[k], x[1]) != 0){ hit_white++; }
					if(cms_[1].at(z[j])(y[k], x[1]) != 0){ hit_grey++; }
				}
				else{ miss ++; }
			}
		}
	}
	if(print){
		std::stringstream ss;
		ss << print_name << "-" << translation[0] << "-" << translation[1] << "-" << translation[2];
		printVectorMatrix(worker, ss.str());
	}

	std::vector<int> r;
	r.push_back(miss);
	r.push_back(hit_white);
	r.push_back(hit_grey);
	r.push_back(tries);
	return r;
}

void Brain::fitMesh(){
	std::vector<std::string> paths;
	paths.push_back("../brain-data/mesh/lh-plial-hull-flood-0-1-merge-4-dof-100k.xml");
	paths.push_back("../brain-data/mesh/rh-plial-hull-flood-0-1-merge-4-dof-100k.xml");
	paths.push_back("../brain-data/mesh/lh-white-hull-flood-0-1-merge-5-dof-100k.xml");
	paths.push_back("../brain-data/mesh/rh-white-hull-flood-0-1-merge-5-dof-100k.xml");

	int eps = 10;
	for(int i = 0; i < paths.size(); i++){
		std::shared_ptr<dolfin::Mesh> mesh = std::make_shared<dolfin::Mesh>(paths.at(i));
		if(rank_ == 0 && verbose_ > 2){	std::cout << "greedy fit mesh " << paths.at(i) << std::endl; }
		auto result = greedyOptimalTranslation(mesh, eps, false, "paths-" + i);
	}

}

std::vector<int> Brain::greedyOptimalTranslation(std::shared_ptr<dolfin::Mesh> mesh, int eps, bool print, std::string print_name){
	if(rank_ == 0 && verbose_ > 2){	std::cout << "find optimal translation from " << print_name <<" to brainweb within eps= " << eps << "..." << std::endl; }
	// update internal mesh
	readMesh(mesh);

	// rough alignment of mesh to brainweb datamap
	int fittedTranslation[3];
	fittedTranslation[0] = min_cm_[0] - (int)floor(min_mesh_[0]);
	fittedTranslation[1] = min_cm_[1] - (int)floor(min_mesh_[1]);
	fittedTranslation[2] = max_cm_[2] - (int)ceil(max_mesh_[2]);
	if(rank_ == 0 && verbose_ > 3){ std::cout << "	bounding translation: " << fittedTranslation[0] << ", " << fittedTranslation[1] << ", " << fittedTranslation[2] << std::endl; }

	// calculate bounds of translation
	int x_start = (((int)floor(min_mesh_[0]) + fittedTranslation[0] - eps) > 0 ? (fittedTranslation[0] - eps) : (int)floor(min_mesh_[0]));
	int y_start = (((int)floor(min_mesh_[1]) + fittedTranslation[1] - eps) > 0 ? (fittedTranslation[1] - eps) : (int)floor(min_mesh_[1]));
	int z_start = (((int)floor(min_mesh_[2]) + fittedTranslation[2] - eps) > 0 ? (fittedTranslation[2] - eps) : (int)floor(min_mesh_[2]));
	int x_end = (((int)ceil(max_mesh_[0]) + fittedTranslation[0] + eps) < dim_[0] ? (fittedTranslation[0] + eps) : (dim_[0] - (int)ceil(max_mesh_[0])));
	int y_end = (((int)ceil(max_mesh_[1]) + fittedTranslation[1] + eps) < dim_[1] ? (fittedTranslation[1] + eps) : (dim_[1] - (int)ceil(max_mesh_[1])));
	int z_end = (((int)ceil(max_mesh_[2]) + fittedTranslation[2] + eps) < dim_[2] ? (fittedTranslation[2] + eps) : (dim_[2] - (int)ceil(max_mesh_[2])));

	// calculate best translation with respect to misses
	std::vector<int> best;
	best.push_back(99999999);			// 	total tries
	best.push_back(99999999);			//	#misses
	best.push_back(-1);					//  x
	best.push_back(-1);					//  y
	best.push_back(-1);					//  z
	// print
	if(rank_ == 0 && verbose_ > 3){ std::cout << "		best translation at missrate " << (int)((double)best.at(1) / (double)best.at(0) * 100) << "% " << std::flush;	}
	for(int x = x_start; x < x_end; x++){
		for(int y = y_start; y < y_end; y++){
			for(int z = z_start; z < z_end; z++){
				int translation[3] = {x, y, z};
				// get result
				auto result = compareTranslated(mesh, translation, false, "shouldnotbe");
				// compare to best
				if(result.at(0) < best.at(1)){
					best.at(0) = result.at(3);
					best.at(1) = result.at(0);
					best.at(2) = x;
					best.at(3) = y;
					best.at(4) = z;
					if(rank_ == 0 && verbose_ > 3){	std::cout << (int)((double)best.at(1) / (double)best.at(0) * 100) << "% " << std::flush; }
				}
			}
		}
	}

	// print to console
	if(rank_ == 0 && verbose_ > 2){
		std::cout << std::endl << "best translation misses " << (int)((double)best.at(1) / (double)best.at(0) * 100) << "%, translation: " <<
				best.at(2) << ", " << best.at(3) << ", " << best.at(4) << std::endl << std::endl;
	}
	// create datafolder with hits to later plot with python
	if(print){
		int translation[3] = {best.at(2), best.at(3), best.at(4)};
		compareTranslated(mesh, translation, true, print_name);
	}
	return best;
}

std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
	Brain::getConcentrationMap(int slice){
	// case when all slices are needed (3D)
	if(slice < 0 || slice > 180){
		std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*> p;
		p = std::make_pair(cms_, dim_);
		return p;
	}
	// case if a specific slice is needed (2D)
	else{
		std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> slices;
		slices.push_back(cms_.at(0).at(slice));
		slices.push_back(cms_.at(1).at(slice));
		std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> single_vec;
		single_vec.push_back(slices);
		std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*> p;
		int* slice_dim = dim_;
		slice_dim[2] = 1;
		p = std::make_pair(single_vec, dim_);
		return p;
	}
}

void Brain::printVectorMatrix(std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> data, std::string name){
	int x = dim_[0];
	int y = dim_[1];
	int z = dim_[2];
	int plane = x*y;
	// create files to convert with python script to picture/ gif
	if(rank_ == 0 && verbose_ > 3){
		std::cout << "	output slice data to: " << dataParentPath_ << "/" << name << std::endl;
	}
	std::string ss = dataParentPath_ + "/" + name;
	mkdir(ss.c_str(), 0777);
	for(int zz = 0; zz < z; zz++){
		// output debug image
		std::stringstream s;
		s << dataParentPath_ << "/" << name << "/slice-" << zz << ".txt";
		std::ofstream out(s.str(), std::ios::trunc);
		for(int p = 0; p < y; p++){
			for(int q = 0; q < x; q++){
				out << (data.at(zz)(p, q) * 255) << ",";
			}
			out << "\n";
		}
		out.close();
	}
}

void Brain::printBooleanSliceForAlphaShape(){
	if(rank_ == 0 && verbose_ > 3){ std::cout << "	create boolean matrix from in BrainWeb tissue data..." << std::endl; }
	std::string greyPath = dataParentPath_ + "/phantom_1.0mm_normal_gry.rawb";
	std::string whitePath = dataParentPath_ + "/phantom_1.0mm_normal_wht.rawb";
	std::string paths[2] = {whitePath, greyPath};
	int x = dim_[0];
	int y = dim_[1];
	int z = dim_[2];
	int plane = x*y;
	char *buf = new char[1];

	// create empty concentration map
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> worker;
	for(int k = 0; k < z; k++){
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> slice = Eigen::MatrixXd::Zero(y, x);
		worker.push_back(slice);
	}
	// collect all points in concentration map
	for(int s = 0; s < 2; s++){
		// read in tissue data
		std::ifstream data(paths[s], std::ios::binary);
		data.read(buf, 1);
		int i = 0;
		int slice_counter = 0;
		int z_curr = 0;
		while(data){
			// information about current voxel
			double x_curr = slice_counter%x;
			double y_curr = slice_counter/x;
			double value = ((double)((unsigned char)*buf));
			// update value map
			if(worker.at(z_curr)(y_curr, x_curr) == 0){
				worker.at(z_curr)(y_curr, x_curr) = (value > 0 ? 1 : 0);
			}
			i++;
			slice_counter++;
			if(slice_counter == plane){
				slice_counter = 0;
				z_curr++;
			}
			data.read(buf, 1);
		}
		data.close();
		if(rank_ == 0 && verbose_ > 3){ std::cout << "		" << i << " data points from " << paths[s] << std::endl; }
	}
	// go through concentration map and output points
	if(rank_ == 0 && verbose_ > 3){
		std::cout << "	output points to csv" << std::endl;
	}
	std::string ss = "../brain-data/brainweb/tissue-location";
	mkdir(ss.c_str(), 0777);
	for(int zz = 0; zz < z; zz++){
		// output slice
		std::stringstream s;
		s << "../brain-data/brainweb/tissue-location/slice-" << zz << ".csv";
		std::ofstream out(s.str(), std::ios::trunc);
		for(int p = 0; p < y; p++){
			for(int q = 0; q < x; q++){
				if(worker.at(zz)(p, q) > 0){ out << q << "," << p << std::endl; }
			}
		}
		out.close();
	}
}
