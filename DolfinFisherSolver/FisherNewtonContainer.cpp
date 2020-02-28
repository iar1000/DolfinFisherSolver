
#include "FisherNewtonContainer.h"
#include "L2Error2D.h"
#include "L2Error3D.h"

FisherNewtonContainer::FisherNewtonContainer(int rank,
		std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> initialCondition,
		std::shared_ptr<dolfin::Expression> D, double rho, double theta, double dt)
{
	rank_ = rank;
	problem_ = std::make_shared<FisherProblem>(rank, mesh, D, rho, dt, theta);
	solver_ = std::make_shared<dolfin::NewtonSolver>();
	mesh_ = mesh;
	initializer_ = initialCondition;
	initialize();
}

FisherNewtonContainer::FisherNewtonContainer(int rank)
{
	if(rank == 0){
		std::cout << "constructor for seperate file initialization not implemented yet" << std::endl;
	}
}

void FisherNewtonContainer::initialize(){
	// stop solver from logging each iteration on it's own
	dolfin::set_log_level(30);
	// set solver parameters
	solver_->parameters["error_on_nonconvergence"] = false; // make sure no error is thrown when not converged
	solver_->parameters["linear_solver"] = "lu";
	solver_->parameters["convergence_criterion"] = "incremental";
	solver_->parameters["maximum_iterations"] = 50;
	solver_->parameters["relative_tolerance"] = 1e-10;
	solver_->parameters["absolute_tolerance"] = 1e-10;
	// get pointers to problem variables
	auto us = problem_->getUs();
	u0_ = us.at(0);				// is going to hold initial fundtion for solver
	u_ = us.at(1);				// is going to hold high precision run result
	u_low_ = us.at(2);			// is going to hold low precision run result
	u_unchanged_ = us.at(3);	// is going to hold the pre-iteration function for possible reset
	dt_ = problem_->getDt();

	// calculate helper parameters for timestep adaption
	p_ = 2;
	if(problem_->getTheta() == 0.5){ p_ = 1; };
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

int FisherNewtonContainer::solve(double t){
	// start tracking new iteration if available
	if(hasTracker_){ tracker_->newIteration(); }

	// solve problem from state u0_ and store in u_
	if(hasTracker_){ tracker_->startTime(); }
	auto results = solver_->solve(*problem_, *u_->vector());
	if(hasTracker_){ tracker_->endTime(); }

	// update u0_ to new state
	*u0_->vector() = *u_->vector();

	// retrieve and save iteration data from solver
	double solverIterations = static_cast<double>(results.first);
	double converged = static_cast<double>(results.second);
	double residual = solver_->residual();
	std::vector<double> data = {t, *dt_, residual, solverIterations, converged};
	if(hasTracker_){ tracker_->addIterationData(data); }

	// finish current iteration
	if(hasTracker_){ tracker_->endIteration(); }

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
	auto results = solver_->solve(*problem_, *u_->vector());
	// retrieve and save low precision iteration data from solver
	double solverIterations = static_cast<double>(results.first);
	double converged = static_cast<double>(results.second);
	double residual = solver_->residual();
	std::vector<double> data = {t, dt, residual, solverIterations, converged};
	if(hasTracker_){ tracker_->addIterationData(data); }
	// save low precision result
	*u_low_->vector() = *u_->vector();

	// solve problem with high precision in two half steps
	*dt_ = dt/2;
	results = solver_->solve(*problem_, *u_->vector());
	*u0_->vector() = *u_->vector();
	results = solver_->solve(*problem_, *u_->vector());

	// u_unchanged is inital function
	// u_low is low low quallity run
	// u is high quality run

	// stop time
	if(hasTracker_){ tracker_->endTime(); }

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
	r = std::make_pair((nabla <= tol), nabla);
	return r;
}

double FisherNewtonContainer::getP(){
	return p_;
}

void FisherNewtonContainer::attachTracker(RuntimeTracker* tracker){
	tracker_ = tracker;
	hasTracker_ = true;
	std::string format = "t, dt, residual, solverIterations, converged";
	tracker_->addIterationFormat(format);
}

void FisherNewtonContainer::output(double t, std::shared_ptr<dolfin::File> pvdFile){
	*pvdFile << std::pair<const dolfin::Function*, double>( u0_.get() , t);
}
