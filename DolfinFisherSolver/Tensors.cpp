
#include "Tensors.h"

// public methods
//////////////////////////
double linear_interpolation(double x, double *values){
	return values[0] * (1.0 - x) + values[1] * x;
}
double bilinear_interpolation(double *x, double *values){
	double worker[2] = {
			linear_interpolation(x[1], &(values[0])),
			linear_interpolation(x[1], &(values[2]))
		};
	return linear_interpolation(x[0], worker);
}
double bilinear_interpolation_v2(double *x, double *values){
	return values[0] * (1 - x[0]) * (1 - x[1]) +
			values[1] * (1 - x[0]) * x[1] +
			values[3] * x[0] * x[1] +
			values[2] * x[0] * (1 - x[1]);

}
double trilinear_interpolation(double *x, double *values){
	double worker[2] = {
			bilinear_interpolation_v2(&(x[0]), &(values[0])),
			bilinear_interpolation_v2(&(x[0]), &(values[4]))
		};
	return linear_interpolation(x[2], worker);
}
double trilinear_interpolation_v2(double *x, double *values){
	double front = values[0] * (1 - x[0]) * (1 - x[1]) +
			values[1] * (1 - x[0]) * x[1] +
			values[3] * x[0] * x[1] +
			values[2] * x[0] * (1 - x[1]);
	double back = values[4] * (1 - x[0]) * (1 - x[1]) +
			values[5] * (1 - x[0]) * x[1] +
			values[7] * x[0] * x[1] +
			values[6] * x[0] * (1 - x[1]);

	return front * (1.0 - x[2]) + back * x[2];
}

// class TensorConstant
/////////////////////////////
TensorConstant::TensorConstant(int rank, double v){
	v_ = v; // value of the tensor
};

std::string TensorConstant::asString(){
	std::stringstream ss;
	ss << "TensorConstant:\n"
				"	value = " << v_ << std::endl;
	return ss.str();
}

void TensorConstant::eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const{
	values[0] = v_;
};


// class TensorSpatial2D
//////////////////////////////
TensorSpatial2D::TensorSpatial2D(int rank, double dcw, double dcg, std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> cms){
	dw_ = dcw;		// Diffusion coefficient for white matter
	dg_ = dcg;		// Diffusion coefficient for grey matter
	cmw_ = cms[0];	// Concentration matrix of white matter
	cmg_ = cms[1];	// Concentration matrix of grey matter
};

std::string TensorSpatial2D::asString(){
	std::stringstream ss;
	ss << "TensorSpatial2D:\n"
				"	D_w = " << dw_ << std::endl <<
				"	D_g = " << dg_ << std::endl <<
				"	coefficient map = vector<Eigen::Matrix>" << std::endl;
	return ss.str();
}

void TensorSpatial2D::eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const{
	double p[2] = {x.data()[0], x.data()[1]};
	double p_frac[2] = {p[0] - floor(p[0]), p[1] - floor(p[1])};

	// calculate value
	double vw[4] = { cmw_(floor(p[0]), floor(p[1])), cmw_(floor(p[0]), ceil(p[1])),
						cmw_(ceil(p[0]), floor(p[1])), cmw_(ceil(p[0]), ceil(p[1])) };
	double vg[4] = { cmg_(floor(p[0]), floor(p[1])), cmg_(floor(p[0]), ceil(p[1])),
							cmg_(ceil(p[0]), floor(p[1])), cmg_(ceil(p[0]), ceil(p[1])) };

	double pw = bilinear_interpolation_v2(p_frac, vw);
	double pg = bilinear_interpolation_v2(p_frac, vg);

	double result = pw * dw_ + pg * dg_;
	values[0] = result;
};


// class TensorSpatial3D
//////////////////////////////
TensorSpatial3D::TensorSpatial3D(int rank, double dcw, double dcg,
		std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> cms, std::vector<int> translation){
	dw_ = dcw;						// Diffusion coefficient for white matter
	dg_ = dcg;						// Diffusion coefficient for grey matter
	translation_ = translation;		// translation vector
	cmw_ = cms[0];					// Vector of concentration matrix of white matter
	cmg_ = cms[1];					// Vector of concentration matrix of grey matter
};

std::string TensorSpatial3D::asString(){
	std::stringstream ss;
	ss << "TensorSpatial3D:\n"
			"	D_w = " << dw_ << std::endl <<
			"	D_g = " << dg_ << std::endl <<
			"	coefficient map = vector<Eigen::Matrix>" << std::endl;
	return ss.str();
}

void TensorSpatial3D::eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const{

	double p[3] = {x.data()[0] + translation_.at(0), x.data()[1] + translation_.at(1), x.data()[2] + translation_.at(2)};
	double p_frac[3] = {p[0] - floor(p[0]), p[1] - floor(p[1]), p[2] - floor(p[2])};

	auto frontW = cmw_.at(floor(p[2]));
	auto backW = cmw_.at(ceil(p[2]));
	double vw[8] = { frontW(floor(p[0]), floor(p[1])), frontW(floor(p[0]), ceil(p[1])),
			frontW(ceil(p[0]), floor(p[1])), frontW(ceil(p[0]), ceil(p[1])),
			backW(floor(p[0]), floor(p[1])), backW(floor(p[0]), ceil(p[1])),
			backW(ceil(p[0]), floor(p[1])), backW(ceil(p[0]), ceil(p[1])) };
	auto frontG = cmg_.at(floor(p[2]));
	auto backG = cmg_.at(ceil(p[2]));
	double vg[8] = { frontG(floor(p[0]), floor(p[1])), frontG(floor(p[0]), ceil(p[1])),
			frontG(ceil(p[0]), floor(p[1])), frontG(ceil(p[0]), ceil(p[1])),
			backG(floor(p[0]), floor(p[1])), backG(floor(p[0]), ceil(p[1])),
			backG(ceil(p[0]), floor(p[1])), backG(ceil(p[0]), ceil(p[1])) };

	double pw = trilinear_interpolation_v2(p_frac, vw);
	double pg = trilinear_interpolation_v2(p_frac, vg);

	values[0] = pw * dw_ + pg * dg_;
};




