#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <dolfin.h>

#include "ProblemSolverContainer.h"
#include "RuntimeTracker.h"
#include "PrintableComponent.h"

// Timestepper for iterating ProblemSolverContainers
class TimeStepper : public PrintableComponent
{
	int rank_;			// MPI rank of instantiating process
	double dt_min_;		// smallest possible timestep
	double dt_max_;		// biggest possible timestep
	double adaptTol_;	// tolerance for timestep adaption via richardson extrapolation
	double adaptSafety_;// safety factorfor timestep adaption via richardson extrapolation

public:
	TimeStepper(int rank, double dt_min, double dt_max, double rTol, double tSafe);
	/*
	 * run desired simulation
	 * type:				1: constant dt
	 * 						2: adaptive dt /in (dt_min, dt_max)
	 * verbosity:			1: no output
	 * 						2: output progress
	 * 						3: output progress and every iteration's information
	 * framesPerTimeUnit:	> 0: number of frames saved to file per time unit
	 * 						= 0: no frames saved as file
	 * 						< 0: all frames saved to file
	*/
	RuntimeTracker run(int type, int verbose,
			std::shared_ptr<dolfin::File> pvdFile, std::string csvPath, int framesPerTimeUnit,
			ProblemSolverContainer* problemContainer, double T, double dt_init);
	std::string asString();	// @override PrintableComponent

private:
	// perform timestepping with constant dt
	// return if convergence failes or duration T is reached
	void constTimestepping(int verbose, double frameDuration, std::shared_ptr<dolfin::File> pvdFile,
			ProblemSolverContainer* problemContainer, double T, double dt_init);
	// perform timestepping with adaptive dt
	// return if T is reached or convergence with dt_min_ still failes
	void adaptiveTimestepping(int verbose, double frameDuration, std::shared_ptr<dolfin::File> pvdFile,
			ProblemSolverContainer* problemContainer, double T, double dt_init);
};

#endif
