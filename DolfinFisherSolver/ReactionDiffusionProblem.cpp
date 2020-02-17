
#include <chrono>

#include "ReactionDiffusionProblem.h"
#include "VariationalReactionDiffusion2D.h"
#include "VariationalReactionDiffusion3D.h"

ReactionDiffusionProblem::ReactionDiffusionProblem(int rank, std::shared_ptr<dolfin::Mesh> mesh, std::shared_ptr<dolfin::Expression> D, double rho, double dt, double theta)
{
	hasTracker_ = false;
	rank_ = rank;
	mesh_ = mesh;
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
	u_temp_ = std::make_shared<dolfin::Function>(V); 	// function holding temporary concentration
	u_ = std::make_shared<dolfin::Function>(V); 	 	// function holding concentration values at t=n+1
	rho_ = std::make_shared<dolfin::Constant>(rho); 	// reaction coefficient
	D_ = D; 											// diffusion tensor
	dt_ = std::make_shared<dolfin::Constant>(dt);		// initial size of time-step
	theta_ = std::make_shared<dolfin::Constant>(theta); // theta of time discretization
	// @TODO: add to ufl function and attach to there
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

void ReactionDiffusionProblem::F(dolfin::GenericVector& b, const dolfin::GenericVector& x){
	auto start = std::chrono::steady_clock::now();
	dolfin::assemble(b, *F_);
	auto end = std::chrono::steady_clock::now();
	if(hasTracker_){
			tracker_->addAssemblyData(std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), 1);
		}
};

void ReactionDiffusionProblem::J(dolfin::GenericMatrix& A, const dolfin::GenericVector& x){
	auto start = std::chrono::steady_clock::now();
	dolfin::assemble(A, *J_);
	auto end = std::chrono::steady_clock::now();
	if(hasTracker_){
		tracker_->addAssemblyData(std::chrono::duration_cast<std::chrono::seconds>(end - start).count(), 2);
	}
};

std::vector<std::shared_ptr<dolfin::Function>> ReactionDiffusionProblem::getUs(){
	return std::vector<std::shared_ptr<dolfin::Function>>({u0_, u_, u_temp_});
}

std::shared_ptr<dolfin::Constant> ReactionDiffusionProblem::getDt(){
	return dt_;
}

std::shared_ptr<dolfin::Mesh> ReactionDiffusionProblem::getMesh(){
	return mesh_;
}

double ReactionDiffusionProblem::getTheta(){
	return *theta_;
}
void ReactionDiffusionProblem::addTracker(RuntimeTracker* tracker, bool record){
	tracker_ = tracker;
	hasTracker_ = record;
	if(rank_ == 0){
		std::cout << "tracker has been added to problem instance, recording = " << record << std::endl;
	}
	// add local dofs to csv
	tracker_->addLocalDofs(V_->dofmap()->index_map()->size(dolfin::IndexMap::MapSize::GLOBAL), V_->dofmap()->index_map()->size(dolfin::IndexMap::MapSize::ALL), rank_);
}

std::string ReactionDiffusionProblem::asString(){
	std::stringstream ss;
	ss << "ReactionDiffusionProblem:" << std::endl <<
					"	mesh = pass by pointer <Mesh> " << std::endl <<
					"	  cell type = " << dolfin::CellType::type2string(mesh_->type().cell_type()) << std::endl <<
					"	  num cells = " << mesh_->num_cells() << std::endl <<
					"	  num vertices = " << mesh_->num_vertices() << std::endl <<
					"	  dof's = " << V_->dofmap()->global_dimension() << std::endl <<
					"	D = pass by pointer <Expression>" << std::endl <<
					"	rho = " << *rho_ << std::endl <<
					"	dt = " << *dt_ << std::endl <<
					"	theta = " << *theta_ << std::endl << std::endl;
	return ss.str();
}


