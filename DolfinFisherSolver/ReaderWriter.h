#ifndef READERWRITER_H
#define READERWRITER_H

#include <dolfin.h>

// divides string at delimiter and returns the substrings as vector
std::vector<std::string> splitString(std::string toSplit, char delim);

// class supporting reading in files/meshes and writing to output
class ReaderWriter
{
	int rank_;								// rank of processor using this ReaderWriter
	std::string outputParent_;				// relative path to parent folder for writing output subfolders and files
	std::string meshParent_;				// relative path to parent folder holding meshes
	std::vector<std::string> components_;	// vector to save strings of component details

public:
	ReaderWriter(int rank, std::string toOutput, std::string toMesh);

	// returns path relativ path to mesh with name "name"
	// returns (file ending/relative path) if success, ("fail"/"fail") otherwise
	std::pair<std::string, std::string> loadMesh(std::string name);

	// creates and returns path for output file with name "name" in outputParent/subfolder
	// if the subdirectory does not exist it is created, must be valid directory name otherwise undefinded behavior
	// type check is done inside function if file type is currently supported
	// returns if creating the path was a success and path if so
	std::pair<bool, std::string> getFilePath(std::string subfolder, std::string filename, std::string suffix);

	// prints the stored list of infos into the subfolder with name _INFO-{filename}
	// thought to be used before a run is started to know what simulation details are
	// use same subfolder as in getFilePath(), otherwise undefinded behavior
	void createRunInfo(std::string subfolder, std::string filename);

	// add string of component details to the list
	void addComponent(std::string);
};

#endif
