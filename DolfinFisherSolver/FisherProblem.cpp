
#include <chrono>

#include "FisherProblem.h"
#include "VariationalFisherEquation2D.h"
#include "VariationalFisherEquation3D.h"

FisherProblem::FisherProblem(int rank, std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> D, double rho, double dt, double theta)
{
	rank_ = rank;
	mesh_ = mesh;
	// create function space, non-linear function F and it's Jacobian J from variational problem definition
	// automatically determine dimensionality of the mesh
	std::shared_ptr<dolfin::FunctionSpace> V;
	std::shared_ptr<dolfin::Form> F, J;
	if(mesh->geometry().dim() == 2){
		V = std::make_shared<VariationalFisherEquation2D::FunctionSpace>(mesh);
		F = std::make_shared<VariationalFisherEquation2D::ResidualForm>(V);
		J = std::make_shared<VariationalFisherEquation2D::JacobianForm>(V, V);
	}
	else if(mesh->geometry().dim() == 3){
		V = std::make_shared<VariationalFisherEquation3D::FunctionSpace>(mesh);
		F = std::make_shared<VariationalFisherEquation3D::ResidualForm>(V);
		J = std::make_shared<VariationalFisherEquation3D::JacobianForm>(V, V);
	}
	else{
		std::cout << "WARNING: only 2 or 3 spatial dimensional meshes allowed" << std::endl;
	}
	// initialize problem parameters
	u0_ = std::make_shared<dolfin::Function>(V); 		// function holding concentration values at t=n
	u_ = std::make_shared<dolfin::Function>(V); 	 	// function holding concentration values at t=n+1
	u_low_ = std::make_shared<dolfin::Function>(V); ;	// function holding concentration values after low resolution run
	u_temp_ = std::make_shared<dolfin::Function>(V); 	// function holding temporary concentration values
	rho_ = std::make_shared<dolfin::Constant>(rho); 	// reaction coefficient
	D_ = D; 											// diffusion tensor
	dt_ = std::make_shared<dolfin::Constant>(dt);		// initial size of time-step
	theta_ = std::make_shared<dolfin::Constant>(theta); // theta of time discretization
	// collect coefficients into groups
	std::map<std::string, std::shared_ptr<const dolfin::GenericFunction>> coefsJ =
		{{"u", u_}, {"rho", rho_}, {"D", D_} , {"dt", dt_}, {"theta", theta_}};
	std::map<std::string, std::shared_ptr<const dolfin::GenericFunction>> coefsF =
			coefsJ;
	coefsF.insert({"u0", u0_});
	// attach coefficients to forms
	J->set_coefficients(coefsJ);
	F->set_coefficients(coefsF);
	// store forms internally
	V_ = V;
	F_ = F; 	// ?
	J_ = J;		// Jacobian of F

	};

void FisherProblem::F(dolfin::GenericVector& b, const dolfin::GenericVector& x){
	dolfin::assemble(b, *F_);
};

void FisherProblem::J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x){
	dolfin::assemble(A, *J_);
};

std::vector<std::shared_ptr<dolfin::Function>> FisherProblem::getUs(){
	return std::vector<std::shared_ptr<dolfin::Function>>({u0_, u_, u_low_, u_temp_});
}

std::shared_ptr<dolfin::Constant> FisherProblem::getDt(){
	return dt_;
}

std::shared_ptr<dolfin::Mesh> FisherProblem::getMesh(){
	return mesh_;
}

double FisherProblem::getTheta(){
	return *theta_;
}

std::string FisherProblem::asString(){
	std::stringstream ss;
	ss << "FisherProblem:" << std::endl <<
					"	mesh = pass by pointer <Mesh> " << std::endl <<
					"	  cell type = " << dolfin::CellType::type2string(mesh_->type().cell_type()) << std::endl <<
					"	  num local cells = " << mesh_->num_cells() << std::endl <<
					"	  num local vertices = " << mesh_->num_vertices() << std::endl <<
					"	  global dof's = " << V_->dofmap()->global_dimension() << std::endl <<
					"	D = pass by pointer <Expression>" << std::endl <<
					"	rho = " << *rho_ << std::endl <<
					"	dt = " << *dt_ << std::endl <<
					"	theta = " << *theta_ << std::endl << std::endl;
	return ss.str();
}


