#ifndef PROBLEMSOLVERCONTAINER_H
#define PROBLEMSOLVERCONTAINER_H

#include "RuntimeTracker.h"

// super-class to use for generic timestepper
class ProblemSolverContainer : public PrintableComponent
{
protected:
	int rank_;
	// tracker variables
	bool hasTracker_ = false;
	RuntimeTracker* tracker_;

public:
	// setup solver with  given parameters: newton and krylov solver
	virtual void initializeSolver(bool verbose, double newtontolrel, double newtontolabs, int newtonmaxiter,
			double krylovtolrel, double krylovtolabs, int krylovmaxiter, std::string ls, std::string pc) = 0;
	// solves the problem and updates internal parameters with results
	// @return:
	virtual int solve(double t, double dt) = 0;
	// solve the problem for time adaptive timestepping
	// via richardson extrapolation with timestep size dt, discretization error toleranze tol
	// @return pair (discretization tolerance met, discretization error nabla)
	virtual std::pair<int, double> solveAdaptive(double t, double dt, double tol) = 0;
	// returns p which depends on the problems discretization parameter theta
	//		theta = 0:   explizit euler discretization, p = 2
	//		theta = 1/2: crank/nicholson discretization, p = 1
	//		theta = 1:	 implicit euler discretization, p = 2
	virtual double getP() = 0;
	// attach tracker object to be able to track iteration data
	// formatting of the csv file must be passed to the tracker
	virtual void attachTracker(RuntimeTracker* tracker) = 0;
	// output the current problem state u0_ to the given file
	// @param: t
	virtual void output(double t, std::shared_ptr<dolfin::File> pvdFile) = 0;
	// output the current function u0 to the file "path+counter.h5"
	virtual void outputFunction(double counter, std::string path) = 0;

};

#endif
