
#include <chrono>

#include "ReactionDiffusionProblem.h"
#include "VariationalReactionDiffusion2D.h"
#include "VariationalReactionDiffusion3D.h"

ReactionDiffusionProblem::ReactionDiffusionProblem(int rank, std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> D, double rho, double dt)
{
	// create function space, non-linear function F and it's Jacobian J from variational problem definition
	// automatically determine dimensionality of the mesh
	std::shared_ptr<dolfin::FunctionSpace> V;
	std::shared_ptr<dolfin::Form> F, J;
	if(mesh->geometry().dim() == 2){
		V = std::make_shared<VariationalReactionDiffusion2D::FunctionSpace>(mesh);
		F = std::make_shared<VariationalReactionDiffusion2D::ResidualForm>(V);
		J = std::make_shared<VariationalReactionDiffusion2D::JacobianForm>(V, V);
	}
	else if(mesh->geometry().dim() == 3){
		V = std::make_shared<VariationalReactionDiffusion3D::FunctionSpace>(mesh);
		F = std::make_shared<VariationalReactionDiffusion3D::ResidualForm>(V);
		J = std::make_shared<VariationalReactionDiffusion3D::JacobianForm>(V, V);
	}
	else{
		std::cout << "WARNING: only 2 or 3 spatial dimensional meshes allowed" << std::endl;
	}
	// initialize problem parameters
	u0_ = std::make_shared<dolfin::Function>(V); 		// function holding concentration values at t=n
	u_ = std::make_shared<dolfin::Function>(V); 	 	// function holding concentration values at t=n+1
	rho_ = std::make_shared<dolfin::Constant>(rho); 	// reaction coefficient
	D_ = D; 											// diffusion tensor
	dt_ = std::make_shared<dolfin::Constant>(dt);		// initial size of time-step
	// collect coefficients into groups
	std::map<std::string, std::shared_ptr<const dolfin::GenericFunction>> coefsJ =
		{{"u", u_}, {"rho", rho_}, {"D", D_} , {"dt", dt_}};
	std::map<std::string, std::shared_ptr<const dolfin::GenericFunction>> coefsF =
			coefsJ;
	coefsF.insert({"u0", u0_});
	// attach coefficients to forms
	J->set_coefficients(coefsJ);
	F->set_coefficients(coefsF);
	// store forms internally
	F_ = F; 	// ?
	J_ = J;		// Jacobian of F
};

void ReactionDiffusionProblem::F(dolfin::GenericVector& b, const dolfin::GenericVector& x){
	auto start = std::chrono::system_clock::now();
	dolfin::assemble(b, *F_);
	auto end = std::chrono::system_clock::now();
	auto duration = end - start;
	std::cout << "assemble F: " << duration.count() / 60 << std::endl;
};

void ReactionDiffusionProblem::J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x){
	auto start = std::chrono::system_clock::now();
	dolfin::assemble(A, *J_);
	auto end = std::chrono::system_clock::now();
	auto duration = end - start;
	std::cout << "assemble J: " << duration.count() / 60 << std::endl;
};

std::pair<std::shared_ptr<dolfin::Function>, std::shared_ptr<dolfin::Function>> ReactionDiffusionProblem::getUs(){
	return std::pair<std::shared_ptr<dolfin::Function>, std::shared_ptr<dolfin::Function>> (u0_, u_);
}

std::string ReactionDiffusionProblem::asString(){
	std::stringstream ss;
	ss << "ReactionDiffusionProblem:" << std::endl <<
					"	mesh = pass by pointer <Mesh> " << std::endl <<
					"	D = pass by pointer <Expression>" << std::endl <<
					"	rho = " << *rho_ << std::endl <<
					"	dt = " << *dt_ << std::endl;
	return ss.str();
}


