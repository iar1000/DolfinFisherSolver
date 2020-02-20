#ifndef SOLVABLEPROBLEM_H
#define SOLVABLEPROBLEM_H

// superclass to solve generic problem-solver couples with the timestepper
// all functions are virtual and must be overwritten by subclass
class SolvableProblem{

	// attaching a tracker to be able to track iteration data
	// what data is tracked exactly is left to define by the instance and has to be told to the tracker
	bool hasTracker = false;
	virtual void attachTracker() = 0;

	// solve the problem from u_{n} to u_{n+1}
	virtual void solve() = 0;

	// solve the problem from u_{n} to u_{n+1} with low precision (one step: dt) and high precision (two steps: 1/2dt)
	// @return: discretization error (nabla), calculated by richardson extrapolation
	virtual double solveRichardson() = 0;

	// @DO THIS INTERNALLY OR IN THE TIMESTEPPER?
	virtual void  outPutToSolutionFunction(double t) = 0;


};

#endif
