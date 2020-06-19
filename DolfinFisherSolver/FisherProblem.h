#ifndef FISHERPROBLEM_H
#define FISHERPROBLEM_H

#include <dolfin.h>
#include "PrintableComponent.h"
#include "RuntimeTracker.h"

class FisherProblem : public dolfin::NonlinearProblem, PrintableComponent
{
	 int rank_;
	 int quad_deg_;								// degree of quadrature
	 std::shared_ptr<dolfin::Mesh> mesh_;		// storing the mesh the function spaces where created from
	 std::shared_ptr<dolfin::FunctionSpace> V_;	// Function space of problem
	 std::shared_ptr<const dolfin::Form> F_;
	 std::shared_ptr<const dolfin::Form> J_;	// Jacobian of F
	 std::shared_ptr<dolfin::Function> u0_;		// function holding concentration values at t=n
	 std::shared_ptr<dolfin::Function> u_;		// function holding concentration values at t=n+1
	 std::shared_ptr<dolfin::Function> u_low_;	// function holding concentration values after low resolution run
	 std::shared_ptr<dolfin::Function> u_temp_;// function holding temporary concentration values
	 std::shared_ptr<dolfin::Constant> rho_;	// reaction coefficient
	 std::shared_ptr<dolfin::Expression> D_; 	// diffusion tensor
	 std::shared_ptr<dolfin::Constant> dt_;		// size of time-step
	 std::shared_ptr<dolfin::Constant> theta_;	// theta value for time discretization

public:
	 FisherProblem(int rank, int quad_deg, std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> D, double rho, double dt, double theta);
	// user defined residual vector
	void F(dolfin::GenericVector& b, const dolfin::GenericVector& x);
	// user defined assamble of Jacobian
	void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x);
	// returns vector of pointers to concentration functions of the problem
	// @return: vector (u0, u, u_low, u_temp)
	std::vector<std::shared_ptr<dolfin::Function>> getUs();
	// returns pointer to dt of variational problem
	std::shared_ptr<dolfin::Constant> getDt();
	// returns theta used for time discretization of problem
	double getTheta();
	// return pointer to fucntion space
	std::shared_ptr<dolfin::FunctionSpace> getV();
	// returns pointer to the mesh used creating function spaces
	std::shared_ptr<dolfin::Mesh> getMesh();
	// @override PrintableComponent
	std::string asString();
};

#endif
