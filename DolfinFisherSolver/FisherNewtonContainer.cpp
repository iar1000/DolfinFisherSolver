
#include "FisherNewtonContainer.h"
#include "L2Error2D.h"
#include "L2Error3D.h"

FisherNewtonContainer::FisherNewtonContainer(int rank,
		std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> initialCondition,
		std::shared_ptr<dolfin::Expression> D, double rho, double theta, double dt)
{
	rank_ = rank;
	problem_ = std::make_shared<FisherProblem>(rank, 4, mesh, D, rho, dt, theta); //
	mesh_ = mesh;
	initializer_ = initialCondition;

	// stop solver from logging each iteration on it's own
	dolfin::set_log_level(40); //16 for trace
	// get pointers to problem variables
	auto us = problem_->getUs();
	u0_ = us.at(0);				// is going to hold initial fundtion for solver
	u_ = us.at(1);				// is going to hold high precision run result
	u_low_ = us.at(2);			// is going to hold low precision run result
	u_unchanged_ = us.at(3);	// is going to hold the pre-iteration function for possible reset
	dt_ = problem_->getDt();

	// calculate helper parameters for timestep adaption
	p_ = 1;
	if(problem_->getTheta() == 0.5){ p_ = 2; };
	Ms_.push_back(L2Error2D::Functional(mesh_, u_low_, u_));
	Ms_.push_back(L2Error3D::Functional(mesh_, u_low_, u_));
	if(mesh_->geometry().dim() == 2){ MIndex_ = 0;}
	else if(mesh_->geometry().dim() == 3){ MIndex_ = 1;}

	// initialize concentration functions
	*u0_ = *initializer_;
	*u_ =  *initializer_;
	*u_unchanged_ = *initializer_;
	*u_low_ = *initializer_;
}

FisherNewtonContainer::FisherNewtonContainer(int rank)
{
	if(rank == 0){
		std::cout << "constructor for seperate file initialization not implemented yet" << std::endl;
	}
}


void FisherNewtonContainer::initializeSolver(bool verbose, double newtontolrel, double newtontolabs, int newtonmaxiter,
		double krylovtolrel, double krylovtolabs, int krylovmaxiter, std::string ls, std::string pc){
	// instanciate krylov solver
	krylovSolver_ = std::make_shared<dolfin::PETScKrylovSolver>(ls, pc);
	newtonSolver_ = std::make_shared<dolfin::NewtonSolver>(mesh_->mpi_comm(), krylovSolver_, dolfin::PETScFactory::instance());
	newtonSolver_->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	newtonSolver_->parameters("krylov_solver")["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	newtonSolver_->parameters["convergence_criterion"] = "residual";
	newtonSolver_->parameters["maximum_iterations"] = newtonmaxiter;
	newtonSolver_->parameters["relative_tolerance"] = newtontolrel;
	newtonSolver_->parameters["absolute_tolerance"] = newtontolabs;
	newtonSolver_->parameters["linear_solver"] = ls;
	newtonSolver_->parameters["preconditioner"] = pc;
	// safe parameters
	ls_ = ls;
	pc_ = pc;
	newton_tol_rel_ = newtontolrel;
	newton_tol_abs_ = newtontolabs;

	if(rank_ == 0 && verbose){	std::cout << newtonSolver_->parameters.str(true) << std::endl; };
}

int FisherNewtonContainer::solve(double t, double dt){
	// start tracking new iteration if available
	if(hasTracker_){ tracker_->newIteration(); }

	// store initial state for possible reset later
	*u_unchanged_->vector() = *u0_->vector();

	// solve problem from state u0_ and store in u_
	*dt_ = dt;
	if(hasTracker_){ tracker_->startTime(); }
	// time the solve kernel
	dolfin::Timer timer("AAA dolfin::NewtonSolver.solve");
	auto results = newtonSolver_->solve(*problem_, *u_->vector());
	timer.stop();
	double wall = std::get<0>(timer.elapsed());	 	// wall time precision around 1 mikro second
	if(hasTracker_){ tracker_->endTime(); }

	// normalize solution
	/*
	double max_u = u_->vector()->max();   				// get the local max
	max_u = dolfin::MPI::max(MPI_COMM_WORLD, max_u);  	// get the global max
	*u_->vector() /= max_u;
	*/

	// retrieve and save iteration data from solver
	double newtonIterations = static_cast<double>(results.first);
	double krylovIterations = static_cast<double>(newtonSolver_->krylov_iterations());
	double converged = static_cast<double>(results.second);
	double residual = newtonSolver_->residual();
	std::vector<double> data = {t, *dt_, residual, newtonIterations, krylovIterations, converged, wall};
	if(hasTracker_){ tracker_->addIterationData(data); }

	// finish current iteration
	if(hasTracker_){ tracker_->endIteration(); }

	// update u0_ depending on convergence
	if(converged){ *u0_->vector() = *u_->vector(); }
	else { *u0_->vector() = *u_unchanged_->vector(); }

	// return convergence info
	return converged;
}


std::pair<int, double> FisherNewtonContainer::solveAdaptive(double t, double dt, double tol){
	// start tracking new iteration if available
	if(hasTracker_){
		tracker_->newIteration();
		tracker_->startTime();
	}

	// store initial state for possible reset later
	*u_unchanged_->vector() = *u0_->vector();

	// solve problem with low precision
	*dt_ = dt;

	// time the solve kernel
	dolfin::Timer timer("AAA dolfin::NewtonSolver.solve");
	auto results = newtonSolver_->solve(*problem_, *u_->vector());
	timer.stop();
	double wall = std::get<0>(timer.elapsed());	 	// wall time precision around 1 mikro second

	// retrieve and save low precision iteration data from solver
	double newtonIterations = static_cast<double>(results.first);
	double krylovIterations = static_cast<double>(newtonSolver_->krylov_iterations());
	double converged_double = static_cast<double>(results.second);
	bool converged = results.second;
	double residual = newtonSolver_->residual();
	std::vector<double> data = {t, dt, residual, newtonIterations, krylovIterations, converged_double, wall};
	if(hasTracker_){ tracker_->addIterationData(data); }

	// sovler converged on the big step, can continue with the small steps
	if(converged){
		// save low precision result
		*u_low_->vector() = *u_->vector();

		// solve problem with high precision in two half steps
		*dt_ = dt/2;
		results = newtonSolver_->solve(*problem_, *u_->vector());
		*u0_->vector() = *u_->vector();
		results = newtonSolver_->solve(*problem_, *u_->vector());

		// stop time
		if(hasTracker_){ tracker_->endTime(); }

		// check convergence
		if(results.second){
			// calculate discretization error nabla
			double errorSqr = dolfin::assemble(Ms_.at(MIndex_));
			double nabla = sqrt(errorSqr) / (pow(2.0,p_) - 1);

			// update state if discretization error tolerance is met
			if(nabla <= tol){ *u0_->vector() = *u_->vector(); }
			// reset state if otherwise
			else{ *u0_->vector() = *u_unchanged_->vector(); }

			// end iteration
			if(hasTracker_){ tracker_->endIteration(); }

			// return discretization criteria decision and nabla
			std::pair<int, double> r;
			r = std::make_pair((nabla <= tol ? 1 : 0), nabla);
			return r;
		}
		else{
			// set nabla such that time step is reduced to 0.5 in timestepper
			double nabla = 2*tol;

			// reset the state
			*u0_->vector() = *u_unchanged_->vector();

			// end iteration
			if(hasTracker_){ tracker_->endIteration(); }

			// return discretization criteria decision and nabla
			std::pair<int, double> r;
			r = std::make_pair(2, nabla);
			return r;
		}
	}
	// big step didn't converge
	else{
		// stop time
		if(hasTracker_){ tracker_->endTime(); }

		// set nabla such that time step is reduced to 0.5 in timestepper
		double nabla = 2*tol;

		// reset the state
		*u0_->vector() = *u_unchanged_->vector();

		// end iteration
		if(hasTracker_){ tracker_->endIteration(); }

		// return discretization criteria decision and nabla
		std::pair<int, double> r;
		r = std::make_pair(2, nabla);
		return r;
	}
}

double FisherNewtonContainer::getP(){
	return p_;
}

std::shared_ptr<FisherProblem> FisherNewtonContainer::getProblem(){
	return problem_;
}

void FisherNewtonContainer::attachTracker(RuntimeTracker* tracker){
	tracker_ = tracker;
	hasTracker_ = true;
	std::string format = "t, dt, residual, newton iterations, krylov iterations, converged, time solve kernel";
	tracker_->addIterationFormat(format);
}

void FisherNewtonContainer::output(double t, std::shared_ptr<dolfin::File> pvdFile){
	*pvdFile << std::pair<const dolfin::Function*, double>( u0_.get() , t);
}

void FisherNewtonContainer::outputFunction(double counter, std::string path){
	std::stringstream ss;
	ss << path << counter << ".h5";
	auto hdf5 = dolfin::HDF5File(MPI_COMM_WORLD, ss.str(), std::string("w"));
	hdf5.write(*u0_, "/u", false);
}

std::string FisherNewtonContainer::asString(){
	std::stringstream ss3;
	ss3 << "FisherNewtonContainer: Solver= " << ls_ << ", Preconditioner= " << pc_ << std::endl <<
			"	Newton tolerance (abs/rel)= (" << newton_tol_abs_ << "/" << newton_tol_rel_ << ")" << std::endl << std::endl <<
			problem_->asString();
	return ss3.str();
}

