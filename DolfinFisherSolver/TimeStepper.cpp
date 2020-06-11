
#include "TimeStepper.h"
#include <sys/stat.h>
#include <fstream>


// class TimeStepper
////////////////////////////
TimeStepper::TimeStepper(int rank, double dt_min, double dt_max, double rTol, double rSafe){
	rank_ = rank;				// MPI rank of instantiating processor
	dt_min_ = dt_min;		// minimum timestep (not used currently)
	dt_max_ = dt_max;		// maximum timestep	(not used currently)
	adaptTol_ = rTol;		// tolerance for timestep adaption via richardson extrapolation
	adaptSafety_ = rSafe;	// safety factorfor timestep adaption via richardson extrapolation
};

std::string TimeStepper::asString(){
	std::stringstream ss;
	ss << "TimeStepper:" << std::endl <<
			"	min dt = " << std::setprecision (15) << dt_min_ << std::endl <<
			"	max dt = " << dt_max_ << std::endl <<
			"	adaptive timestep tolerance = " << adaptTol_ << std::endl <<
			"	adaptive timestep factor = " << adaptSafety_ << std::endl;
	return ss.str();
}

RuntimeTracker TimeStepper::run(int simulationType, int verbose,
		std::shared_ptr<dolfin::File> pvdFile, std::string csvPath, int framesPerTimeUnit,
		ProblemSolverContainer* problemContainer, double T, double dt_init, int dt_runlength)
{
	// calculate time between two frames
	// frameDuration is multiplied with #frames and compared to t decide when to take next frame
	double frameDuration;
	if(framesPerTimeUnit > 0){ frameDuration = 1.0 / framesPerTimeUnit;	}
	else if(framesPerTimeUnit == 0){ frameDuration = std::numeric_limits<double>::max(); }
	else { frameDuration = 0.0; }

	// instantiate runtime information instance
	// only create rank 0 process tracker with actual verbosity, rest is silent
	RuntimeTracker tracker;
	std::stringstream ss;
	ss << csvPath << "-" << rank_;
	if(rank_ == 0){	tracker = RuntimeTracker(verbose, true, ss.str()); }
	else { tracker = RuntimeTracker(0, false, ss.str());	}

	// create base function file path from csv-path
	// can easily create a bug at some point
	std::string s = ss.str();
	std::string delim = "/";
	auto pos = s.find(delim);
	auto pre_pos = s.find(delim);
	while (pos != std::string::npos)
	{
		pre_pos = pos;
		pos = s.find(delim, pos + 1);
	}
	std::stringstream folderPath;
	folderPath << s.substr(0, pre_pos) << "/functions";
	if(mkdir(folderPath.str().c_str(), 0777) == 0){
		if(rank_ == 0){
			std::cout << "INFO (getFilePath): created subfolder " << folderPath.str() << std::endl;
		}
	}

	std::stringstream functionPath;
	functionPath << s.substr(0, pre_pos) << "/functions/function";


	// attach tracker to problem container, which handles the rest
	problemContainer->attachTracker(&tracker);

	// print runtime information, verbose level 2
	if(rank_ == 0 && verbose > 1){
		std::cout << "run-info: " << std::endl <<
		"	simulation type: " << simulationType << std::endl <<
		"	T = " << T << std::endl <<
		"	dt init = " << dt_init << std::endl <<
		"	csv path = " << csvPath << std::endl <<
		"	function output path = " << functionPath.str() << std::endl;}

	// start tracking simulation
	tracker.startTracking();

	// run correct simulation type
	if(simulationType == 1){ constTimestepping(verbose, frameDuration, pvdFile, functionPath.str(), problemContainer, T, dt_init); }
	else if(simulationType == 2){ adaptiveTimestepping(verbose, frameDuration, pvdFile, functionPath.str(), problemContainer, T, dt_init, dt_runlength); }
	else{ if(rank_ == 0){ std::cout << "WARNING (Timestepper): unknown simulation type, no run!" << std::endl; }}

	// stop tracking simulation
	tracker.stopTracking();

	// return tracker of the simulation
	return tracker;
}

void TimeStepper::constTimestepping(int verbose, double frameDuration, std::shared_ptr<dolfin::File> pvdFile, std::string functionOutPath,
 		ProblemSolverContainer* problemContainer, double T, double dt_init)
{
	// set timestepping and helper variables
	double t = 0;					// current time
	double dt = dt_init;			// current timestep size
	int frameNumber = 0;			// counter for outputed frames
	int functionNumber = 0;			// counter for outputed functions
	int progress = 0;				// progress in %
	double onePercent = T/100; 	// duration of one % of the whole progress

	// print start of timestepping
	if(rank_ == 0){ std::cout << std::endl << "run temporal constant timestepping..." << std::endl << std::endl; }

	// perform timestepping
	while(t < T){

		// run iteration, continue on successful convergence
		int converged = problemContainer->solve(t, dt);
		if(!converged){
			if(rank_ == 0){ std::cout << std::endl << "iteration failed to converge at t = " << t << std::endl; }
			break;
		}

		// update t
		t += dt;

		// write frame
		if(t >= frameNumber * frameDuration){
			// print verbose level 3
			if(rank_ == 0 && verbose > 2){ std::cout << "writing system state to pvd..." << std::endl << std::endl; }
			problemContainer->output(t, pvdFile);
			frameNumber++;
		}

		// write function
		if(t >= functionNumber){
			if(rank_ == 0 && verbose > 2){ std::cout << "writing system state to function..." << std::endl << std::endl; }
			problemContainer->outputFunction(functionNumber, functionOutPath);
			functionNumber++;
		}

		// print progress verbose level 2
		if(rank_ == 0 && verbose > 1 && t >= progress * onePercent){
			std::cout << "simulation status: [" << progress << "/100]%" << std::endl << std::endl;
			progress++;
		}
	}
	// always save last state
	problemContainer->output(t, pvdFile);
	problemContainer->outputFunction(0, functionOutPath);

	// print termination of timestepping
	if(rank_ == 0){ std::cout << "finished temporal constant timestepping!" << std::endl << std::endl; }
}

void TimeStepper::adaptiveTimestepping(int verbose, double frameDuration, std::shared_ptr<dolfin::File> pvdFile, std::string functionOutPath,
		ProblemSolverContainer* problemContainer, double T, double dt_init, int dt_runlength)
{
	// set timestepping and helper variables
	double t = 0;						// current time
	double dt = dt_init;				// current timestep size
	double dtNew = dt_init;				// adapted timestep size
	int frameNumber = 0;				// counter for outputed frames
	int functionNumber = 0;				// counter for outputed functions
	int progress = 0;					// progress in %
	double onePercent = T/100;			// duration of one % of the whole progress
	double p = problemContainer->getP();// get theta dependend p for timestep adaption


	// print start of timestepping
	if(rank_ == 0){ std::cout << std::endl << "run temporal adaptive timestepping with runlength " << dt_runlength << "..." <<  std::endl << std::endl; }

	// perform timestepping
	problemContainer->output(t, pvdFile);
	while(t < T){

		// run iteration
		auto results = problemContainer->solveAdaptive(t, dt, adaptTol_);
		bool criteriumMet = results.first;
		double nabla = results.second;

		// update timestep
		double fac = adaptSafety_ * pow((adaptTol_ / nabla), (1/p));

		// cap the maximum increase
		if(fac > 1.01){
			fac = 1.01;
		}

		dtNew = fac * dt;

		// finish iteration if discretization error tolerance is met
		// state in ProblemSolverContainer should be updated
		if(criteriumMet){
			// print discretation tolerance met, verbose level 3
			if(rank_ == 0 && verbose > 2){
				std::cout << "time discretiation error criterium met\n	increase dt to "<< dtNew << "  (" << fac << ")" << std::endl << std::endl;
			}

			// update t
			t += dt;

			// write frame
			if(t >= frameNumber * frameDuration){
				// print verbose level 3
				if(rank_ == 0 && verbose > 2){ std::cout << "writing system state to pvd..." << std::endl << std::endl; }
				problemContainer->output(t, pvdFile);
				frameNumber++;
			}

			// write function
			if(t >= functionNumber){
				if(rank_ == 0 && verbose > 2){ std::cout << "writing system state to function..." << std::endl << std::endl; }
				problemContainer->outputFunction(functionNumber, functionOutPath);
				functionNumber++;
			}

			// print progress verbose level 2
			if(rank_ == 0 && verbose > 1 && t >= progress * onePercent){
				std::cout << "simulation status: [" << progress << "/100]%" << std::endl;
				progress++;
			}
		}
		// don't update t, discretization error tolerance was not met
		// problemContainer didn't update internal state, so next run will solve same system with adapted dt
		else{
			// print discretization tolerance not met and details, verbose level 3
			if(rank_ == 0 && verbose > 2){
				std::cout << "time discretiation error tolerance not met: " << std::endl <<
						"	discretization error (nabla) = " << nabla << " ( > " << adaptTol_ << ")" << std::endl <<
						"	decrease of dt: " << dt << " -> " << dtNew << " (" << (1.0/fac) << ")" << std::endl << std::endl;
			}
		}
		// update dt
		dt = dtNew;

	}
	// always save last state
	problemContainer->output(t, pvdFile);
	problemContainer->outputFunction(0, functionOutPath);

	// print termination of timestepping
	if(rank_ == 0){ std::cout << "finished temporal adaptive timestepping!" << std::endl << std::endl; }
}





