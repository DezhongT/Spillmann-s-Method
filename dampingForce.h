#ifndef DAMPINGFORCE_H
#define DAMPINGFORCE_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"

class dampingForce
{
public:
	dampingForce(elasticRod &m_rod, timeStepper &m_stepper, double m_viscosity, double m_eta_per, double m_eta_par);
	~dampingForce();
	void computeFd();
	void computeJd();
	void computeFd1();

private:
	elasticRod *rod;
	timeStepper *stepper;
	double viscosity, eta_per, eta_par;
    double headForceCoeff;
	double headSize;
	
    Vector3d t, u, f;
    int ind, indx, indy;
    Matrix3d Id3, jac;
    double dt;
    Matrix<double,1,3> v;
};

#endif
