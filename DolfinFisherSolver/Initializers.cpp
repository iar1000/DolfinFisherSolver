
#include "Initializers.h"
#include <math.h>

// class InitializeCircle
/////////////////////////////
InitializerCircle::InitializerCircle(double x, double y, double r, double v){
		x0_ = x; // x-coordinate
		y0_ = y; // y-coordinate
		r0_ = r; // radius
		v0_ = v; // initializer value
};
std::string InitializerCircle::asString(){
	std::stringstream ss;
	ss << "InitializerCircle:" << std::endl <<
			"	radius = " << r0_ << std::endl <<
			"	value = " << v0_ << std::endl <<
			"	location = (" << x0_ << "," << y0_ << ")" << std::endl;
	return ss.str();
}
void InitializerCircle::eval(dolfin::Array<double> &values, const dolfin::Array<double> &x) const
{
    if ((x[0] - x0_) * (x[0] - x0_) + (x[1] - y0_) * (x[1] - y0_) <= r0_ * r0_){
    	values[0] = v0_;
    }
    else { values[0] = 0; }
};


// class InitializerSphere
////////////////////////////////
InitializerSphere::InitializerSphere(double x, double y, double z, double r, double v){
		x0_ = x; // x-coordinate
		y0_ = y; // y-coordinate
		z0_ = z; // z-coordinate
		r0_ = r; // radius
		v0_ = v; // initializer value
};
std::string InitializerSphere::asString(){
	std::stringstream ss;
	ss << "InitializerSphere:" << std::endl <<
			"	radius = " << r0_ << std::endl <<
			"	value = " << v0_ << std::endl <<
			"	location = (" << x0_ << "," << y0_ << "," << z0_ << ")" << std::endl;
	return ss.str();
}
void InitializerSphere::eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const
{
	if ((x[0] - x0_) * (x[0] - x0_) + (x[1] - y0_) * (x[1] - y0_) + (x[2] - z0_) * (x[2] - z0_) <= r0_ * r0_){
		values[0] = v0_;
    }
	else { values[0] = 0; }
};


// class InitializerGaussianSphere
////////////////////////////////
InitializerGaussianSphere::InitializerGaussianSphere(double x, double y, double z, double r, double v){
		x0_ = x; // x-coordinate
		y0_ = y; // y-coordinate
		z0_ = z; // z-coordinate
		r0_ = r; // radius of normal distribution
		v0_ = v; // maximum cell density
};
std::string InitializerGaussianSphere::asString(){
	std::stringstream ss;
	ss << "InitializerGaussianSphere:" << std::endl <<
			"	radius of normal distribution= " << r0_ << std::endl <<
			"	maximum value = " << v0_ << std::endl <<
			"	location = (" << x0_ << "," << y0_ << "," << z0_ << ")" << std::endl;
	return ss.str();
}
void InitializerGaussianSphere::eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const
{
	double dist = (x[0] - x0_) * (x[0] - x0_) + (x[1] - y0_) * (x[1] - y0_) + (x[2] - z0_) * (x[2] - z0_);
	double value_raw = exp(-dist / r0_);
	values[0] = v0_ * value_raw;
};

