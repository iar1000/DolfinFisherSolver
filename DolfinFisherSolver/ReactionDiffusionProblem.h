#ifndef IOHandler_H
#define IOHandler_H

#include <dolfin.h>
#include "PrintableComponent.h"

class ReactionDiffusionProblem : public dolfin::NonlinearProblem, PrintableComponent
{
	 std::shared_ptr<const dolfin::Form> F_;
	 std::shared_ptr<const dolfin::Form> J_;	// Jacobian of F
	 std::shared_ptr<dolfin::Function> u0_;		// function holding concentration values at t=n
	 std::shared_ptr<dolfin::Function> u_;		// function holding concentration values at t=n+1
	 std::shared_ptr<dolfin::Constant> rho_;	// reaction coefficient
	 std::shared_ptr<dolfin::Expression> D_; 	// diffusion tensor
	 std::shared_ptr<dolfin::Constant> dt_;		// initial size of time-step

public:
	ReactionDiffusionProblem(int rank, std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> D, double rho, double dt);
	// user defined residual vector
	void F(dolfin::GenericVector& b, const dolfin::GenericVector& x);
	// user defined assamble of Jacobian
	void J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x);
	// returns pair of pointers to concentration functions
	// @return: pair (u0, u)
	std::pair<std::shared_ptr<dolfin::Function>, std::shared_ptr<dolfin::Function>> getUs();
	// @override PrintableComponent
	std::string asString();
};

#endif
