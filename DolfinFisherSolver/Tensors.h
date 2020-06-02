#ifndef TENSORS_H
#define TENSORS_H

#include <vector>
#include <dolfin.h>
#include <Eigen/Dense>

#include "PrintableComponent.h"

// linear interpolation
// interpolates at target point x [0, 1] between v0 and v1
double linear_interpolation(double x, double *values);

// bilinear interpolation
// interpolates target point (x[0], x[1])
// values given are in (x,y) order 0,0 - 0,1 - 1,0 - 1,1
double bilinear_interpolation(double *x, double *values);

// direct calculated linear interpolation
// interpolates target point (x[0], x[1])
// values given are in (x,y) order 0,0 - 0,1 - 1,0 - 1,1
double bilinear_interpolation_v2(double *x, double *values);

// trilinear interpolation
// interpolates target point (x[0], x[1], x[2])
// values are given in (x,y,z) order defining two faces
//				front [0-3]: 0,0,0 - 0,1,0 - 1,0,0 - 1,1,0
//				back  [4-7]: 0,0,1 - 0,1,1 - 1,0,1 - 1,1,1
double trilinear_interpolation(double *x, double *values);

// trilinear interpolation, without recursion
// interpolates target point (x[0], x[1], x[2])
// values are given in (x,y,z) order defining two faces
//				front [0-3]: 0,0,0 - 0,1,0 - 1,0,0 - 1,1,0
//				back  [4-7]: 0,0,1 - 0,1,1 - 1,0,1 - 1,1,1
double trilinear_interpolation_v2(double *x, double *values);

// Hashable coordinate-entry for interpolation buffer
struct HashableCoordinates
{
	double x;
	double y;
	double z;
	HashableCoordinates(double xx, double xy, double xz){
		x = xx; y = xy; z = xz;
	}
	bool operator==(const HashableCoordinates &other) const{
		return (x == other.x
				&& y == other.y
				&& z == other.z);
	}
	friend std::ostream& operator<<(std::ostream& os, const HashableCoordinates& h)
	{
	    os << "[" << h.x << ',' << h.y << ',' << h.z << "]";
	    return os;
	}
};
struct CoordinateHasher
{
	size_t operator()( const HashableCoordinates& k ) const
	{
		// Compute individual hash values for first, second and third
		// http://stackoverflow.com/a/1646913/126995
		std::size_t res = 17;
		res = res * 31 + std::hash<double>()( k.x );
		res = res * 31 + std::hash<double>()( k.y );
		res = res * 31 + std::hash<double>()( k.z );
		return res;
	}
};

// Constant Tensor
// Expression returning same value for all evaluated cells
class TensorConstant : public dolfin::Expression, PrintableComponent
{
	double v_;	// value of the tensor
	int rank_;
public:
	std::string asString(); 	// @override PrintableComponent
	TensorConstant(int rank, double v);
	void eval (dolfin::Array<double> &values, const dolfin::Array<double> &x, const ufc::cell &cell) const override;
};


// 2D spatial dependent Tensor
// Expression returning value dependent on the spatial location of the cell
class TensorSpatial2D : public dolfin::Expression, PrintableComponent
{
	int rank_;
	double dw_;		// Diffusion coefficient for white matter
	double dg_;		// Diffusion coefficient for grey matter
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cmw_;	// Concentration matrix of white matter
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cmg_;	// Concentration matrix of grey matter
	std::unordered_map<HashableCoordinates, double, CoordinateHasher> interpolationBuffer_; // stores already interpolated values

public:
	std::string asString(); 	// @override PrintableComponent
	TensorSpatial2D(int rank, double dcw, double dcg, std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> cms);
	void eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const override;
	void addBufferEntry(HashableCoordinates p, double res);
};


// 3D spatial dependent Tensor
// Expression returning value dependent on the spatial location of the cell
class TensorSpatial3D : public dolfin::Expression, PrintableComponent
{
	bool useBuffer;	// For performance comparison between buffer and non-buffer
	double dw_;		// Diffusion coefficient for white matter
	double dg_;		// Diffusion coefficient for grey matter
	std::vector<int> translation_; // translation size x, y, z which is added to mesh coordinates to fit on brainweb matrices
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> cmw_;	// Vector of concentration matrix of white matter
	std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> cmg_;	// Vector of concentration matrix of grey matter
	std::unordered_map<HashableCoordinates, double, CoordinateHasher> interpolationBuffer_; // stores already interpolated values
public:
	std::string asString(); 	// @override PrintableComponent
	TensorSpatial3D(int rank, double dcw, double dcg,
			std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> cms,
			std::vector<int> translation);
	void eval (dolfin::Array<double> &values, const dolfin::Array<double> &x, const ufc::cell &cell) const override;
	void addBufferEntry(HashableCoordinates p, double res);
};


#endif
