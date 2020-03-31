
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <sys/stat.h>
#include <sstream>


#include "Brain.h"

Brain::Brain(int verbose, std::string dataParentPath){
	verbose_ = verbose;
	dataParentPath_ = dataParentPath;
	readData();
}

void Brain::readData(){
	if(verbose_ > 3){ std::cout << "read in BrainWeb tissue data..." << std::endl; }
	std::string greyPath = dataParentPath_ + "/phantom_1.0mm_normal_gry.rawb";
	std::string whitePath = dataParentPath_ + "/phantom_1.0mm_normal_wht.rawb";
	double range = 255;
	int x = 181;
	int y = 217;
	int z = 181;
	dim_[0] = x;
	dim_[1] = y;
	dim_[2] = z;
	int plane = x*y;
	char *buf = new char[1];

	// read in white tissue
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> worker1;
	for(int k = 0; k < z; k++){
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> slice = Eigen::MatrixXd::Zero(y, x);
		worker1.push_back(slice);
	}
	std::ifstream white_data(whitePath, std::ios::binary);
	white_data.read(buf, 1);
	int i = 0;
	int slice_counter = 0;
	int slice = 0;
	while(white_data){
		worker1.at(slice)(slice_counter/x, slice_counter%x) = ((double)((unsigned char)*buf)) / range;
		i++;
		slice_counter++;
		if(slice_counter == plane){
			slice_counter = 0;
			slice++;
		}
		white_data.read(buf, 1);
	}
	white_data.close();
	cms_.push_back(worker1);
	if(verbose_ > 3){ std::cout << "	done reading " << i << " white tissue data points" << std::endl; }

	// read in grey tissue
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> worker2;
	for(int k = 0; k < z; k++){
		Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> slice = Eigen::MatrixXd::Zero(y, x);
		worker2.push_back(slice);
	}
	std::ifstream grey_data(greyPath, std::ios::binary);
	grey_data.read(buf, 1);
	i = 0;
	slice_counter = 0;
	slice = 0;
	while(grey_data){
		worker2.at(slice)(slice_counter/x, slice_counter%x) = ((double)((unsigned char)*buf)) / range;
		i++;
		slice_counter++;
		if(slice_counter == plane){
			slice_counter = 0;
			slice++;
		}
		grey_data.read(buf, 1);
	}
	grey_data.close();
	cms_.push_back(worker2);
	if(verbose_ > 3){ std::cout << "	done reading " << i << " grey tissue data points" << std::endl; }

	// create files to convert with python script to picture/ gif
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



}

std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*>
	Brain::getConcentrationMap(){
	std::pair<std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>, int*> p;
	p = std::make_pair(cms_, dim_);
	return p;
}
