
#include <sys/stat.h>
#include <fstream>

#include "ReaderWriter.h"

// public methods
//////////////////////////////////////
std::vector<std::string> splitString(std::string toSplit, char delim){
	std::vector<std::string> tokens;
	std::stringstream ss(toSplit);
	std::string token;
	while(std::getline(ss, token, delim)){
		tokens.push_back(token);
	}
	return tokens;
}

// class ReaderWriter
//////////////////////////////////
ReaderWriter::ReaderWriter(int rank, std::string toOutput, std::string toMesh){
	rank_ = rank;
	outputParent_ = toOutput;
	meshParent_ = toMesh;

	if(rank_ == 0){
		std::cout << "created ReaderWriter:" << std::endl <<
				"	path to output parent = " << outputParent_ << std::endl <<
				"	path to mesh parent = " << meshParent_ << std::endl;
	}
}
std::pair<bool, int> ReaderWriter::loadMesh(std::shared_ptr<dolfin::Mesh> mesh, std::string name)
{
	// check input format
	std::vector<std::string> tokens = splitString(name, '.');
	if(tokens.size() != 2){
		if(rank_ == 0){
			std::cout << "WARNING: wrong mesh name format" << std::endl;
		}
		return std::pair<bool,int>(false, 0);
	}
	// check if source format of mesh is supported
	std::string format = tokens.at(1);
	if(format == "h5"){
		try{
			auto hdf5 = dolfin::HDF5File(MPI_COMM_WORLD, meshParent_ + "/" + name, std::string("r"));
			hdf5.read(*mesh, "/mesh", false);
			// add mesh to components list before returning
			std::stringstream ss;
			ss << "Mesh:" << std::endl <<
					"	name = " << name << std::endl <<
					"	dimension = " << mesh->geometry().dim() << std::endl;
			addComponent(ss.str());

			return std::pair<bool,int>(true, mesh->geometry().dim());
		}
		catch(const std::exception& e){
			if(rank_ == 0){
				std::cout << "WARNING: mesh not found. path: " << meshParent_ << "/" << name << std::endl;
			}
			return std::pair<bool,int>(false, 0);
		}
	}
	else{
		if(rank_ == 0){
			std::cout << "WARNING: couldn't load mesh: " << name <<
					". datatype not supported (yet)" << std::endl;
		}
		return std::pair<bool,int>(false, 0);
	}
}

std::pair<bool, std::string> ReaderWriter::getFilePath(std::string subfolder, std::string filename, std::string suffix){
	// check if output file type is supportet
	std::string folderPath = outputParent_ + "/" + subfolder;
	if(suffix == "pvd" ||		// output format for paraview
			suffix == "txt"		// output format for infofile
					){
		// create sub-directory
		if(mkdir(folderPath.c_str(), 0777) == 0){
			if(rank_ == 0){
				std::cout << "INFO (getFilePath): created subfolder " << folderPath << std::endl;
			}
		}
		else{
			if(errno == EEXIST){
				if(rank_ == 0){
				std::cout << "INFO (getFilePath): already existing subfolder " << folderPath << std::endl;
				}
			}
			else{
				if(rank_ == 0){
					std::cout << "INFO (getFilePath): failed creating subfolder " << folderPath << std::endl;
				}
				return std::pair<bool, std::string>(false, folderPath + "/" + filename + "." + suffix);
			}
		}

		// return path
		std::string path = folderPath + "/" + filename + "." + suffix;
		return std::pair<bool, std::string>(true, path);
	}
	// file format no supported
	else{
		if(rank_ == 0){
			std::cout << "WARNING: no support for output format (yet) " << suffix << std::endl;
		}
		return std::pair<bool, std::string>(false, folderPath + "/" + filename + "." + suffix);
	}
}

void ReaderWriter::createRunInfo(std::string subfolder, std::string filename){
	// create path to infofile
	auto pathReturn = getFilePath(subfolder, "_INFO-" + filename, "txt");
	if(pathReturn.first){
		std::string path = pathReturn.second;
		std::ofstream file(path);
		// print components
		for ( auto i = components_.begin(); i != components_.end(); i++ ) {
		    file << *i << std::endl;
		}
		file.close();
	}
	else{
		std::cout << "WARNING: no infofile created at " << pathReturn.second << std::endl;
	}
}

void ReaderWriter::addComponent(std::string details){
	components_.push_back(details);
}
