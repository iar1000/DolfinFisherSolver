#ifndef INITIALIZERS_H
#define INITIALIZERS_H

#include <dolfin.h>

#include "PrintableComponent.h"

// Initializer used for 2D meshes
// Expression to set the value of all cells in the input circle to a given value
class InitializerCircle : public dolfin::Expression, PrintableComponent
{
	double x0_; // x-coordinate
	double y0_;	// y-coordinate
	double r0_; // radius
	double v0_; // initializer value
public:
	std::string asString();	//@override PrintableCOmponent
	InitializerCircle(double x, double y, double r, double v);
	void eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const override;
};


// Initializer used for 3D meshes
// Expression to set the value of all cells in the input sphere to a given value
class InitializerSphere : public dolfin::Expression, PrintableComponent
{
	double x0_; // x-coordinate
	double y0_;	// y-coordinate
	double z0_;	// z-coordinate
	double r0_; // radius
	double v0_; // initializer value
public:
	std::string asString();	//@override PrintableCOmponent
	InitializerSphere(double x, double y, double z, double r, double v);
	void eval (dolfin::Array<double> &values, const dolfin::Array<double> &x) const override;
};

#endif
