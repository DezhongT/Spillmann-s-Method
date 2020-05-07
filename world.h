#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"

// include external force
#include "dampingForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

//include collision checker
#include "collision.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void OpenFile1(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	void updateTimeStep_data();

	void CoutDataC(ofstream &outfile);
		
private:

	// Physical parameters
	double RodLength;
	double helixradius, helixpitch;
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;

	// double EI;
	// double GJ;

	// double youngM;

	double pulltime;
	double friction;
	double v_constant;
	double statictime;


	// Viscous drag coefficients
    double eta_per, eta_par;
    
	double tol, stol;
	int maxIter; // maximum number of iterations
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	VectorXd theta;
	
	// Rod
	elasticRod *rod;
	
	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;
	double currentTime;
	
	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	elasticTwistingForce *m_twistingForce;
	inertialForce *m_inertialForce;	
	externalGravityForce *m_gravityForce;
	dampingForce *m_dampingForce;

	collision *check;

	
	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();

	void updateBoundary();

	void updateCons();

	void newtonMethod(bool &solved);

	void newtonMethodC(bool &solved);
	void calculateForce();
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	bool hold;


	vector<int> hold_for_tying;


	vector<int> contactZ;

	Vector3d temp;
	Vector3d temp1;
	Vector3d gravity;
	Vector3d inertial;
	Vector3d dampingF;

	int iter_c;
};

#endif
