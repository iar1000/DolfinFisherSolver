#ifndef IOHandler_H
#define IOHandler_H

#include <dolfin.h>
#include "PrintableComponent.h"
#include "TimeStepper.h"

class RuntimeTracker; 	// forward definition of runtime tracker b.c. of cyclic dependency

class ReactionDiffusionProblem : public dolfin::NonlinearProblem, PrintableComponent
{
	 std::shared_ptr<dolfin::Mesh> mesh_;		// storing the mesh the function spaces where created from
	 std::shared_ptr<dolfin::FunctionSpace> V_;	// Function space of problem
	 std::shared_ptr<const dolfin::Form> F_;
	 std::shared_ptr<const dolfin::Form> J_;	// Jacobian of F
	 std::shared_ptr<dolfin::Function> u0_;		// function holding concentration values at t=n
	 std::shared_ptr<dolfin::Function> u_;		// function holding concentration values at t=n+1
	 std::shared_ptr<dolfin::Function> u_temp_; // function for temporarily storing
	 std::shared_ptr<dolfin::Constant> rho_;	// reaction coefficient
	 std::shared_ptr<dolfin::Expression> D_; 	// diffusion tensor
	 std::shared_ptr<dolfin::Constant> dt_;		// size of time-step
	 std::shared_ptr<dolfin::Constant> theta_;	// theta value for time discretization
	 int rank_;
	 bool hasTracker_;							// bool indicating if tracker is initialized
	 RuntimeTracker* tracker_;					// tracker for recording time data


public:
	ReactionDiffusionProblem(int rank, std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> D, double rho, double dt, double theta);
	// user defined residual vector
	void F(dolfin::GenericVector& b, const dolfin::GenericVector& x);
	// user defined assamble of Jacobian
	void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x);
	// returns vector of pointers to concentration functions
	// @return: vector (u0, u_mid, u)
	std::vector<std::shared_ptr<dolfin::Function>> getUs();
	// returns pointer to dt of variational problem
	std::shared_ptr<dolfin::Constant> getDt();
	// returns theta used for time discretization of problem
	double getTheta();
	// returns pointer to the mesh used creating function spaces
	std::shared_ptr<dolfin::Mesh> getMesh();
	// adds tracker and option to start recording data
	void addTracker(RuntimeTracker* tracker, bool record);
	// @override PrintableComponent
	std::string asString();
};

#endif
