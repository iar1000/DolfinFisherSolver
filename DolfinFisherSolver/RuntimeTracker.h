
#ifndef RUNTIMETRACKER_H
#define RUNTIMETRACKER_H

#include <dolfin.h>
#include <chrono>
#include <vector>

#include <fstream>

#include "PrintableComponent.h"

// class for tracking all the information of a run
class RuntimeTracker : public PrintableComponent
{
	int verbose_;								// verbosity of runtime tracker
	int numberIterations_;						// number successfully calculated iterations
	std::vector<std::vector<double>> iterations_;									// buffer for tracked iterations
	std::chrono::time_point<std::chrono::system_clock> startTimeRun_, endTimeRun_;		// start and stop time of simulation
	int elapsedAll_;							// elapsed time in seconds between start and stop of simulation
	// iteration state
	std::string format_;						// format of iteration data delemited by ","
	bool inIteration_;							// if an iteration is currently tracked
	std::vector<double> currIteration_;			// If inIteration, stores added iteration data
	std::chrono::time_point<std::chrono::system_clock> startTimeIteration_;	// start time of simulation, end time is filled automatically into currIteration
	int elapsedIt_;								// time elapsed in the current iteration
	// storing iteration data
	bool toCsv_;								// if this tracker outputs its data to csv
	int iterationBufferSize_ = 50;			// number of iteration stored before writing to file
	std::string csvPath_;						// path to file storing iteration data as csv
	std::ofstream csv_;							// csv file for iteration data

public:
	std::string asString();	// @override PrintableComponent
	RuntimeTracker(int verbose, bool toCsv, std::string csvPath);
	RuntimeTracker();
	void addIterationFormat(std::string format); // sets format for iteration data, use "," delemiter between elements
	// tracking simulation
	void startTracking();	// sets start time of whole simulation
	void stopTracking();	// sets end time of whole simulation
	// tracking iteration
	void newIteration();	// prepare tracker for a fresh iteration
	void startTime();		// take time when started solver
	void endTime();			// take time when finished solver
	void addIterationData(std::vector<double> data); // enter data for current iteration
	void endIteration();	// saves all data collected between last newIteration() and this function call
};

#endif
