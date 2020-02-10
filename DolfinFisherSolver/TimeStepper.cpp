
#include "TimeStepper.h"
#include "L2Error2D.h"


// class RuntimeInformation
///////////////////////////////
RuntimeTracker::RuntimeTracker(){}
RuntimeTracker::RuntimeTracker(int simulationType, int verbose, double T, double dt_min, double dt_max){
	// general state
	verbose_ = verbose;					// verbosity of tracker
	simulationType_ = simulationType;	// type of simulation running
	T_ = T;								// simulation time
	dt_min_ = dt_min;					// minimum timestep
	dt_max_ = dt_max;					// maximum timestep
	elapsedAll_ = -1;					// if elapsedAll_ > 0, tracking has been ended
	// iterations state
	inIteration_ = false;				// no running iteration currently
	numberIterations_ = 0;				// number of tracked iterations

}
std::string RuntimeTracker::asString(){
	std::stringstream ss;
	// decide type
	if(simulationType_ == 1){
		ss << "RuntimeTracker:" << std::endl <<
				"	simulation type = constant time" << std::endl <<
				"	initializer = pass by pointer <Expression>" << std::endl <<
				"	output = pass by pointer <File>" << std::endl;
	}
	else if(simulationType_ == 2){
		ss << "RuntimeTracker:" << std::endl <<
				"	simulation type = adaptive time" << std::endl <<
				"	initializer = pass by pointer <Expression>" << std::endl <<
				"	output = pass by pointer <File>" << std::endl;
	}
	// tracking has been ended, add simulation recap
	if(elapsedAll_ > 0){
		int minutes = elapsedAll_ / 60;
		ss << std::endl << "	sucessfully run " << (iterations_.back().t / T_ * 100) << " % ("<< numberIterations_ << " iterations)" << std::endl <<
				"	elapsed time: " << (minutes / 60) << "h " << (minutes % 60) << " min" << std::endl;
	}
	return ss.str();
}

void RuntimeTracker::startTracking(){
	startTimeRun_ = std::chrono::system_clock::now();
	// print
	if(verbose_ > 1) {
		std::time_t now = std::chrono::system_clock::to_time_t(startTimeRun_);
		std::string time(30, '\0');
		std::strftime(&time[0], time.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
		std::cout << "start tracker: (" << time << ")" << std::endl;
	}
}
void RuntimeTracker::stopTracking(){
	endTimeRun_ = std::chrono::system_clock::now();
	std::chrono::duration<double> duration = endTimeRun_ - startTimeRun_;
	elapsedAll_ = duration.count();
	// print
	if(verbose_ > 1) { // verbose level 2
		std::time_t now = std::chrono::system_clock::to_time_t(endTimeRun_);
		std::string time(30, '\0');
		std::strftime(&time[0], time.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
		std::cout << "stop tracker: (" << time << ")" << std::endl <<
				"	elapsed time: " << elapsedAll_/60 << " min" << std::endl <<
				"	elapsed time: " << elapsedAll_/3600 << " h" << std::endl;
	}
}
// tracking iteration
void RuntimeTracker::newIteration(){
	inIteration_ = true;
	currIteration_ = Iteration();
}
void RuntimeTracker::endIteration(){
	if(inIteration_){
		inIteration_ = false;
		// iterations_.push_back(currIteration_); @TODO: do i need to track this data?
		numberIterations_++;
		// print
		if(verbose_ > 2){ // verbose level 3
			std::cout << "	t = " << currIteration_.t << " (total iteration " << numberIterations_ << ")" << std::endl <<
					"	size of dt = " << currIteration_.dt << std::endl <<
					"	converged = " << (currIteration_.converged ? "true" : "false") << std::endl <<
					"	residual = " << currIteration_.residual << std::endl <<
					"	solver iterations = " << currIteration_.numberSteps << std::endl <<
					"	elapsed time = " << currIteration_.elapsedTime << std::endl << std::endl;
		}
	}
	else{
		std::cout << "WARNING: seems like no iteration has been started yet!" << std::endl <<
				"	No information stored" << std::endl;
	}
}
void RuntimeTracker::startTime(){
	startTimeIteration_ = std::chrono::system_clock::now();
}
void RuntimeTracker::endTime(){
	auto endTimeIteration = std::chrono::system_clock::now();
	std::chrono::duration<double> duration = endTimeIteration - startTimeIteration_;
	currIteration_.elapsedTime = duration.count();
}
void RuntimeTracker::addIterationData(double t, double dt, bool converged, int numSteps, double residual){
	currIteration_.t = t;
	currIteration_.dt = dt;
	currIteration_.converged = converged;
	currIteration_.numberSteps = numSteps;
	currIteration_.residual = residual;
}


// class TimeStepper
////////////////////////////
TimeStepper::TimeStepper(int rank,
		std::shared_ptr<ReactionDiffusionProblem> problem,	std::shared_ptr<dolfin::NewtonSolver> solver,
		double T, double dt_min, double dt_max){
	rank_ = rank;				// MPI rank of instantiating processor
	problem_ = problem;			// Nonlinear problem to solve
	solver_ = solver; 			// newton solver
	solver_->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	dolfin::set_log_level(40);								// stop solver from logging
	T_ = T;					// simulation time
	dt_min_ = dt_min;		// minimum timestep
	dt_max_ = dt_max;		// maximum timestep
};

std::string TimeStepper::asString(){
	std::stringstream ss;
	ss << "TimeStepper:" << std::endl <<
			"	problem = pass by pointer <ReactionDiffusionProblem>" << std::endl <<
			"	solver = pass by pointer <NewtonSolver>" << std::endl <<
			"	T = " << T_ << std::endl <<
			"	min dt = " << std::setprecision (15) << dt_min_ << std::endl <<
			"	max dt = " << dt_max_ << std::endl;
	return ss.str();
}

RuntimeTracker TimeStepper::run(int simulationType, int verbose, std::shared_ptr<dolfin::Expression> initializer,
		std::shared_ptr<dolfin::File> output, int framesPerTimeUnit, double dt)
{
	// decide frame duration
	// frameDuration is multiplied with #frames to decide when to take next frame, if left 0 all frames are taken
	double frameDuration = 0.0;
	if(framesPerTimeUnit > 0){
		frameDuration = 1.0 / framesPerTimeUnit;
	}

	// instantiate runtime information instance
	// only create rank 0 process tracker with actual verbosity, rest is silent
	RuntimeTracker tracker;
	if(rank_ == 0){
		tracker = RuntimeTracker(simulationType, verbose, T_, dt_min_, dt_max_);
	}
	else {
		tracker = RuntimeTracker(simulationType, 0, T_, dt_min_, dt_max_);
	}
	// print runtime information
	if(rank_ == 0){
		std::cout << tracker.asString();
	}

	tracker.startTracking();	// start tracker

	// decide type and run simulation
	if(simulationType == 1){
		constTime_timestepping(verbose, &tracker, initializer, output, frameDuration, dt);
	}
	else if(simulationType == 2){
		adaptiveTime_timestepping(verbose, &tracker, initializer, output, frameDuration, dt);
	}
	else{
		if(rank_ == 0){
			std::cout << "WARNING: unknown simulation type" << std::endl;
		}
	}

	tracker.stopTracking();

	// return runtime information of simulation
	return tracker;
}

void TimeStepper::constTime_timestepping(int verbose, RuntimeTracker *tracker, std::shared_ptr<dolfin::Expression> initializer, std::shared_ptr<dolfin::File> output, double frameDuration, double tdt)
{
	// set timestepping and helper variables
	double t = 0;
	double dt = tdt;
	int frameNumber = 0;
	// get and initialize concentration function
	auto us = problem_->getUs();
	std::shared_ptr<dolfin::Function> u0 = us.at(0);
	std::shared_ptr<dolfin::Function> u = us.at(1);
	*u0 = *initializer;
	*u =  *initializer;

	if(rank_ == 0){ std::cout << "running..." << std::endl << std::endl; }
	int progress = 0;			// progress in %
	double onePercent = T_/100; // one % of the whole progress

	*output << std::pair<const dolfin::Function*, double>( u.get() , t); 	// save initial state of the system
	while(t < T_){
		// prepare current iteration
		tracker->newIteration();
		t += dt;
		*u0->vector() = *u->vector();

		// solve problem and track results
		tracker->startTime();
		auto results = solver_->solve(*problem_, *u->vector());
		tracker->endTime();

		int solverIterations = results.first;
		bool converged = results.second;
		double residual = solver_->residual();
		//tracker->addIterationData(t, dt, converged, solverIterations, residual);

		// check convergence
		if(!converged){
			// finalize this iteration and quit
			tracker->endIteration();
			break;
		}

		// write frame
		if(t > frameNumber * frameDuration){
			*output << std::pair<const dolfin::Function*, double>( u.get() , t);
			frameNumber++;

			if(rank_ == 0 && verbose > 2){ // verbose level 3
				std::cout << "writing system state to output..." << std::endl << std::endl;
			}
		}

		// print
		if(rank_ == 0 && verbose > 1 && t >= progress * onePercent){  // verbose level 2
			std::cout << "simulation status: [" << progress << "/100]%" << std::endl;
			progress++;
		}


		// finalize this iteration
		tracker->endIteration();
	}
	*output << std::pair<const dolfin::Function*, double>( u.get() , t);  	// save final state of the system

	if(rank_ == 0){ std::cout << "finished" << std::endl << std::endl; }
}

void TimeStepper::adaptiveTime_timestepping(int verbose, RuntimeTracker *tracker, std::shared_ptr<dolfin::Expression> initializer, std::shared_ptr<dolfin::File> output, double frameDuration, double tdt)
{
	/*
	 * u_n (V)
	 * u_n_low(V)
	 * u_n_high(V)
	 *
	 * solver_low(u_n, u_n_low)
	 * solver_low -> t, dt
	 *
	 * solver_high_1(u_n, u_n_high)
	 * solver_high_1 -> t, dt/2
	 *
	 * solver_high_2(u_n_high, u_n_high)
	 * solver_high_2 -> t+dt, dt/2
	 *
	 * est = compute_est(u_n_low, u_n_high)
	 * dt_new = ...
	 */

	// set timestepping and helper variables
	double t = 0;
	double dt = tdt;		// current timestep
	double dtNew = tdt;		// next timestep size
	int frameNumber = 0;
	// get pointers to problem variables
	auto us = problem_->getUs();
	std::shared_ptr<dolfin::Function> u0_p = us.at(0);
	std::shared_ptr<dolfin::Function> u_p = us.at(1);
	std::shared_ptr<dolfin::Function> u_low = us.at(2);
	auto dt_p = problem_->getDt();
	// initialize L2-Norm functional and helpers for timestep adaption
	double tol = 10e-8;
	double safety = 0.9;
	std::shared_ptr<dolfin::Mesh> mesh = problem_->getMesh();
	L2Error2D::Functional M(mesh, u_low, u_p);
	double p = 1;
	if(problem_->getTheta() == 0.5){ p = 2; };
	// initialize concentration function
	*u0_p = *initializer; 		// problem u^n concentration
	*u_p =  *initializer;		// problem u^{n+1} concentration
	*u_low = *initializer;	// temporary concentration

	if(rank_ == 0){ std::cout << "running..." << std::endl << std::endl; }
	int progress = 0;			// progress in %
	double onePercent = T_/100; // one % of the whole progress

	*output << std::pair<const dolfin::Function*, double>( u_p.get() , t); 	// save initial state of the system
	while(t < T_){
		tracker->newIteration();

		// solve problems and track time/results
		tracker->startTime();

		// prepare and run low precision run
		*dt_p = dt;
		auto results = solver_->solve(*problem_, *u_p->vector());
		int solverIterations = results.first;
		bool converged = results.second;
		double residual = solver_->residual();
		tracker->addIterationData(t, dt, converged, solverIterations, residual);
		// u_temp_p holds low end condition
		*u_low->vector() = *u_p->vector();

		// prepare and run high precision run
		*dt_p = dt/2;
		// compute first half of timestep
		results = solver_->solve(*problem_, *u_p->vector());
		// update initial condition to halftime result
		*u0_p->vector() = *u_p->vector();
		// compute second half of timestep
		results = solver_->solve(*problem_, *u_p->vector());
		// u_p holds high end condition

		// calculate richardson extrapolation nabla
		double errorSqr = dolfin::assemble(M);
		double error = sqrt(errorSqr);
		double nabla = error / (pow(2.0,p) - 1);
		double fac = pow((safety * tol / nabla), (1/p));
		dtNew = fac * dt;
		if(nabla > tol){
			if(rank_ == 0 && verbose > 2){	// verbose level 3
				std::cout << "new timestep= " << dtNew << " (" << dt << ")" << std::endl;
			}
			dt = dtNew;	// update to smaller timestep
			continue; // don't update t, repeat with smaller timestep
		}

		// prepare next iteration
		t += dt;
		dt = dtNew;	// update to bigger timestep
		*u0_p->vector() = *u_p->vector(); //@Which solver iteration values are used?

		tracker->endTime();


		// check convergence
		if(!converged){
			// finalize this iteration and quit
			tracker->endIteration();
			break;
		}

		// write frame
		if(t > frameNumber * frameDuration){
			*output << std::pair<const dolfin::Function*, double>( u_p.get() , t);
			frameNumber++;

			if(rank_ == 0 && verbose > 2){ // verbose level 3
				std::cout <<"writing system state to output..." << std::endl << std::endl;
			}
		}

		// print
		if(rank_ == 0 && verbose > 1 && t >= progress * onePercent){  // verbose level 2
			std::cout << "simulation status: [" << progress << "/100]%" << std::endl;
			progress++;
		}

		// finalize this iteration
		tracker->endIteration();
	}
	*output << std::pair<const dolfin::Function*, double>( u_p.get() , t);  	// save final state of the system

	if(rank_ == 0){ std::cout << "finished" << std::endl << std::endl; }

}





