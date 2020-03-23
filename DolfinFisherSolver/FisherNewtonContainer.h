#ifndef FISHERNEWTONCONTAINER_H
#define FISHERNEWTONCONTAINER_H

#include <dolfin.h>

#include "ProblemSolverContainer.h"
#include "FisherProblem.h"

class FisherNewtonContainer : public ProblemSolverContainer
{
	// problem specific variables
	std::shared_ptr<dolfin::Expression> initializer_;	// initial condition of the concentration
	std::shared_ptr<FisherProblem> problem_;			// fisher problem class
	std::shared_ptr<dolfin::Function> u0_;				// pointer to concentration function of fisher-problem t=n
	std::shared_ptr<dolfin::Function> u_;				// pointer to concentration function of fisher-problem t=n+1
	std::shared_ptr<dolfin::Function> u_unchanged_;		// pointer to placeholder storing concentration function at beginning of a adaptive run
	std::shared_ptr<dolfin::Function> u_low_;			// pointer to placeholder for low precision result
	std::shared_ptr<dolfin::Constant> dt_;				// pointer to constant holding size of time-step of fisher-problem
	std::shared_ptr<dolfin::Mesh> mesh_;				// pointer to mesh with which the function space of the problem was created
	int meshDimensions_;									// mesh dimensionality
	double p_;											// p for calculating discretization error (nabla) and new timestep
	std::vector<dolfin::Form> Ms_;						// initialize all dimension forms in the beginning
	int MIndex_;										// decide later which to use
	// solver specific variables
	std::shared_ptr<dolfin::NewtonSolver> newtonSolver_;
	std::shared_ptr<dolfin::PETScKrylovSolver> krylovSolver_;

	// mesh convergence study
	double last_t;
	double sample_range;

public:
	// constructor with parameters for fisher equation
	FisherNewtonContainer(int rank,
			std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> initialCondition,
			std::shared_ptr<dolfin::Expression> D, double rho, double theta, double dt);
	// constructor loading parameters from seperate file
	FisherNewtonContainer(int rank);

	// @override virtuals from ProblemSolverContainer
	void initializeSolver(bool verbose, double newtontolrel, double newtontolabs, int newtonmaxiter,
				double krylovtolrel, double krylovtolabs, int krylovmaxiter, std::string ls, std::string pc);
	int solve(double t);
	std::pair<int, double> solveAdaptive(double t, double dt, double tol);
	double getP();
	void attachTracker(RuntimeTracker* tracker);
	void output(double t, std::shared_ptr<dolfin::File> pvdFile);

};

#endif
