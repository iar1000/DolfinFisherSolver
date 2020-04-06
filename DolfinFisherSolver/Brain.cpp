
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
	for(int i = 0; i < 3; i++){	max_cm_[i] = -99999; min_cm_[i] = 999999;	}
	readData();
}



void Brain::readData(){
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

	// create files to convert with python script to picture/ gif
	/*
	if(verbose_ > 3){
		std::cout << "	output slice data to: " << dataParentPath_ << "/brain_slices_csv" << std::endl;
		std::string ss = dataParentPath_ + "/brain_slices_csv";
		mkdir(ss.c_str(), 0777);
		for(int zz = 0; zz < z; zz++){
			// output debug image
			std::stringstream s;
			s << dataParentPath_ << "/brain_slices_csv/white-z-" << zz << ".txt";
			std::ofstream out(s.str(), std::ios::trunc);
			for(int p = 0; p < y; p++){
				for(int q = 0; q < x; q++){
					out << (worker1.at(zz)(p, q) * 255) << ",";
				}
				out << "\n";
			}
			out.close();
		}
		for(int zz = 0; zz < z; zz++){
			// output debug image
			std::stringstream s;
			s << dataParentPath_ << "/brain_slices_csv/grey-z-" << zz << ".txt";
			std::ofstream out(s.str(), std::ios::trunc);
			for(int p = 0; p < y; p++){
				for(int q = 0; q < x; q++){
					out << (worker2.at(zz)(p, q) * 255) << ",";
				}
				out << "\n";
			}
			out.close();
		}
	}
	*/


}

void Brain::fitMesh(std::shared_ptr<dolfin::Mesh> mesh){
	// variables
	double max[3];
	double min[3];
	for(int i = 0; i < 3; i++){
		max[i] = -99999;
		min[i] = 99999;
	}
	// go through coordinates
	int num_dim = 3;
	auto coordinates = mesh->coordinates();
	for(int i = 0; i < coordinates.size(); i += num_dim){
		// update extremas
		if(coordinates[i] > max[0]){ max[0] = coordinates[i]; }
		if(coordinates[i] < min[0]){ min[0] = coordinates[i]; }
		if(coordinates[i+1] > max[1]){ max[1] = coordinates[i+1]; }
		if(coordinates[i+1] < min[1]){ min[1] = coordinates[i+1]; }
		if(coordinates[i+2] > max[2]){ max[2] = coordinates[i+2]; }
		if(coordinates[i+2] < min[2]){ min[2] = coordinates[i+2]; }
	}
	// extract important parameters
	double lengths[3];
	for(int i = 0; i < 3; i++){ lengths[i] = max[i] - min[i]; }

	if(rank_ == 0 && verbose_ > 3){
		std::cout << "	mesh bounds: " << std::endl <<
				"		x: [" << min[0] << ", " << max[0] << "]" << " (" << lengths[0] << ")" << std::endl <<
				"		y: [" << min[1] << ", " << max[1] << "]" << " (" << lengths[1] << ")" << std::endl <<
				"		z: [" << min[2] << ", " << max[2] << "]" << " (" << lengths[2] << ")" << std::endl;
	}

}

std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
	Brain::getConcentrationMap(){
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*> p;
	p = std::make_pair(cms_, dim_);
	return p;
}
