
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
}

std::pair<std::string, std::string> ReaderWriter::loadMesh(std::string name)
{
	// check input format
	std::vector<std::string> tokens = splitString(name, '.');
	if(tokens.size() != 2)
	{
		if(rank_ == 0){ std::cout << "WARNING (ReaderWriter): wrong mesh name format" << std::endl; }
		return std::pair<std::string,std::string>("fail", "fail");
	}
	// check if source format of mesh is supported
	std::string format = tokens.at(1);
	if(format == "h5"){ return std::pair<std::string,std::string>("h5", meshParent_ + "/" + name); }
	else if(format == "xml"){ return std::pair<std::string,std::string>("xml", meshParent_ + "/" + name); }
	else if(format == "xdmf"){ return std::pair<std::string,std::string>("xdmf", meshParent_ + "/" + name); }
	else { return std::pair<std::string,std::string>("fail", "fail"); }
}

std::pair<bool, std::string> ReaderWriter::getFilePath(std::string subfolder, std::string filename, std::string suffix){
	// check if output file type is supportet
	if(	suffix == "txt"	||		// output format for infofile
			suffix == "csv"		// output format for iteration data
					){
		// create sub-directory
		std::string folderPath = outputParent_ + "/" + subfolder;
		if(mkdir(folderPath.c_str(), 0777) == 0){
			if(rank_ == 0){
				std::cout << "INFO (getFilePath): created subfolder " << folderPath << std::endl;
			}
		}
		else{
			if(errno == EEXIST){

			}
			else{
				if(rank_ == 0){
					char cwd[PATH_MAX];
					std::cout << "Current working dir: " << getcwd(cwd, sizeof(cwd)) << std::endl;
					std::cout << "INFO (getFilePath): failed creating subfolder " << folderPath << std::endl;
				}
				return std::pair<bool, std::string>(false, folderPath + "/" + filename + "." + suffix);
			}
		}

		// return path
		std::string path = folderPath + "/" + filename + "." + suffix;
		return std::pair<bool, std::string>(true, path);
	}
	else if(suffix == "pvd"){		// output format for paraview, create folder b.c. of pvd spam
		// create sub-directory with pvd folder
		std::string folderPath = outputParent_ + "/" + subfolder;
		if(mkdir(folderPath.c_str(), 0777) == 0){
			folderPath = outputParent_ + "/" + subfolder + "/pvd";
			if(mkdir(folderPath.c_str(), 0777) == 0){
				if(rank_ == 0){ std::cout << "INFO (getFilePath): created subfolder " << folderPath << std::endl; }
			}
		}
		else{
			if(errno == EEXIST){
				folderPath = outputParent_ + "/" + subfolder + "/pvd";
				if(mkdir(folderPath.c_str(), 0777) == 0){
					if(rank_ == 0){ std::cout << "INFO (getFilePath): created subfolder " << folderPath << std::endl; }
				}
			}
			else{
				if(rank_ == 0){
					std::cout << "INFO (getFilePath): failed creating subfolder " << folderPath << std::endl;
					std::cout << "ERRNO: " << errno << std::endl;
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
		return std::pair<bool, std::string>(false, "false");
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

void ReaderWriter::updateStatusFile(int status, std::string subfolder, int nprocs){
	std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::string time(30, '\0');
	std::strftime(&time[0], time.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
	// create path to statusfile
	auto pathReturn = getFilePath(subfolder, "_STATUS", "txt");
	if(pathReturn.first){
		std::string path = pathReturn.second;
		// check the status and act accordingly
		if(status == 0){
			std::ofstream file;
			file.open(path, std::ios_base::out | std::ios_base::trunc);
			file << "number processors: " << nprocs << std::endl <<
					"Simulation started - " << time << std::endl;
			file.close();
		}
		else if(status == 1){
			std::ofstream file;
			file.open(path, std::ios_base::app);
			file << "Mesh loaded - " << time << std::endl;
			file.close();
		}
		else if(status == 2){
			std::ofstream file;
			file.open(path, std::ios_base::app);
			file << "Simulation running - " << time << std::endl;
			file.close();
		}
		else if(status == 3){
			std::ofstream file;
			file.open(path, std::ios_base::app);
			file << "Simulation finnished - " << time << std::endl;
			file.close();
		}
		else{
			std::cout << "WARNING: no statusfile created, status not known: " << status << std::endl;
		}
	}
	else{
		std::cout << "WARNING: no statusfile created, path error: " << pathReturn.second << std::endl;
	}
}


void ReaderWriter::addComponent(std::string details){
	components_.push_back(details);
}
