#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <dolfin.h>
#include <vector>
#include <chrono>
#include <fstream>

#include "ReactionDiffusionProblem.h"
#include "PrintableComponent.h"

// class for tracking all the information of a run
class RuntimeTracker : public PrintableComponent
{
	// nested class for storing iteration data
	struct Iteration{
		int elapsedTime;	// elapsed time between startTime() and endTime()
		bool converged;		// if solver converged this iteration
		int numberSteps;	// number of steps of solver till convergence
		double residual;	// residual of solution
		double t;			// simulated timestep of simulation
		double dt;			// last size of timestep used before convergence
	};
	int verbose_;								// verbosity of runtime tracker
	int simulationType_;						// type of simulation running
	double T_;									// simulation time
	double dt_min_;								// minimum timestep
	double dt_max_;								// maximum timestep
	int numberIterations_;						// number successfully calculated iterations
	std::vector<RuntimeTracker::Iteration> iterations_;	// vector storing information about single iterations
	std::chrono::time_point<std::chrono::system_clock> startTimeRun_, endTimeRun_;		// start and stop time of simulation
	int elapsedAll_;							// elapsed time in seconds between start and stop
	// iteration state
	bool inIteration_;							// if an iteration is currently tracked
	RuntimeTracker::Iteration currIteration_;	// If inIteration, Iterationobject to be filled and added to list after finnishing iteration
	std::chrono::time_point<std::chrono::system_clock> startTimeIteration_;
	// storing iteration data
	bool toCsv_;								// if this tracker outputs its data to csv
	int iterationBufferSize_ = 1000;			// number of iteration stored before writing to file
	std::string csvPath_;						// path to file storing iteration data as csv
	std::ofstream csv_;							// csv file for iteration data

public:
	std::string asString();	// @override PrintableComponent
	RuntimeTracker(int simulationType, int verbose, bool toCsv, double T, double dt_min, double dt_max, std::string outputPath);
	RuntimeTracker();
	// tracking simulation
	void startTracking();	// sets start time of whole simulation
	void stopTracking();	// sets end time of whole simulation
	// tracking iteration
	void newIteration();	// prepare tracker for a fresh iteration
	void startTime();		// take time when started solver
	void endTime();			// take time when finished solver
	void addIterationData(double t, double dt, bool converged, int numSteps, double residual); // enter data for current iteration
	void endIteration();	// saves all data collected between last newIteration() and this function call
};

// Timestepper for ReactionDiffusionProblem/ NewtonSolver combination, handling the time-stepping
class TimeStepper : public PrintableComponent
{
	int rank_;			// MPI rank of instantiating process
	std::shared_ptr<ReactionDiffusionProblem> problem_; // Reactiondiffusion problem to solve
	std::shared_ptr<dolfin::NewtonSolver> solver_;		// newton solver
	double T_;			// simulation time
	double dt_min_;		// smallest possible timestep
	double dt_max_;		// biggest possible timestep
	double richTol_;	// tolerance for richardson extrapolation adaptive timestep
	double richSafety_;	// safety factor for richardson extrapolation adaptive timestep

public:
	std::string asString();	// @override PrintableComponent
	TimeStepper(int rank,
			std::shared_ptr<ReactionDiffusionProblem> problem, std::shared_ptr<dolfin::NewtonSolver> solver,
			double T, double dt_min, double dt_max, double rTol, double tSafe);
	/*
	 * run desired simulation
	 * type:				1: constant dt
	 * 						2: adaptive dt /in (dt_min, dt_max)
	 * verbosity:			1: no output
	 * 						2: output progress
	 * 						3: output progress and every iteration's information
	 * framesPerTimeUnit:	>  0:	number of frames save as file per time unit
	 * 						<= 0:	all frames saved as file
	*/
	RuntimeTracker run(int type, int verbose, std::shared_ptr<dolfin::Expression> initializer,
			std::shared_ptr<dolfin::File> output, std::string csvPath,
			int framesPerTimeUnit, double dt);

private:
	// perform timestepping with constant dt
	// return if convergence failes or duration T is reached
	void constTime_timestepping(int verbose, RuntimeTracker *tracker, std::shared_ptr<dolfin::Expression> initializer, std::shared_ptr<dolfin::File> output, double frameDuration, double tdt);
	// perform timestepping with adaptive dt
	// return if T is reached or convergence with dt_min_ still failes
	void adaptiveTime_timestepping(int verbose, RuntimeTracker *tracker, std::shared_ptr<dolfin::Expression> initializer, std::shared_ptr<dolfin::File> output, double frameDuration, double tdt);

};

#endif
