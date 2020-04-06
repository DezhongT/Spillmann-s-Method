#include "world.h"
#include <sstream>
#include <iomanip>

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				// boolean
	saveData = m_inputData.GetBoolOpt("saveData");			// boolean

	// Physical parameters
	RodLength = m_inputData.GetScalarOpt("RodLength");      // meter
    helixradius = m_inputData.GetScalarOpt("helixradius");  // meter
    gVector = m_inputData.GetVecOpt("gVector");             // m/s^2
    maxIter = m_inputData.GetIntOpt("maxIter");             // maximum number of iterations
    helixpitch = m_inputData.GetScalarOpt("helixpitch");    // meter
	rodRadius = m_inputData.GetScalarOpt("rodRadius");      // meter
	numVertices = m_inputData.GetIntOpt("numVertices");     // int_num
	youngM = m_inputData.GetScalarOpt("youngM");            // Pa
	Poisson = m_inputData.GetScalarOpt("Poisson");          // dimensionless
	deltaTime = m_inputData.GetScalarOpt("deltaTime");      // seconds
	totalTime= m_inputData.GetScalarOpt("totalTime");       // seconds
	tol = m_inputData.GetScalarOpt("tol");                  // small number like 10e-7
	stol = m_inputData.GetScalarOpt("stol");				// small number, e.g. 0.1%
	density = m_inputData.GetScalarOpt("density");          // kg/m^3
	viscosity = m_inputData.GetScalarOpt("viscosity");      // viscosity in Pa-s
	pulltime = m_inputData.GetScalarOpt("pulltime");        //get time of pulling
    statictime = m_inputData.GetScalarOpt("statictime");
    friction = m_inputData.GetScalarOpt("friction");
 
	shearM = youngM/(2.0*(1.0+Poisson));					// shear modulus
	
	// Viscous drag coefficients using Resistive Force Theory
	eta_per = 4.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) + 0.5);
    eta_par = 2.0*M_PI*viscosity/( log(2*helixpitch/rodRadius) - 0.5);
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;
    ostringstream num;
    num << std::setprecision(2);
    num << std::fixed;
    num << friction;
    name << "datafiles/simDER" <<num.str()<<".txt";
	outfile.open(name.str().c_str());
	outfile.precision(10);	
	// outfile << "# x [meter] y [meter] z [meter]\n";
}


void world::OpenFile1(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	time_t current_time = time(0);

	// Open an input file named after the current time
	ostringstream name;
    ostringstream num;
    num << std::setprecision(2);
    num << std::fixed;
    num << friction;
    name << "datafiles/simDER_v" <<num.str()<<".txt";
	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
		return;

	outfile.close();
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
		return;

	for (int i = 0; i < rod->nv; i++)
	{
		if (i<rod->ne)
		{
			outfile  << rod->x(4*i) << " " <<rod->
		x(4*i+1) << " " << rod->x(4*i+2) <<" "<< rod->x(4*i+3)<<endl;
		}
		else
		{
			outfile<< rod->x(4*i) << " " <<rod->
		x(4*i+1) << " " << rod->x(4*i+2) <<" "<< 0<<endl;
		}
	}
}


void world::CoutDataC(ofstream &outfile)
{
	double f;
	double f1;

    f = temp.norm();
    f1 = temp1.norm();

	double e;
	temp = rod->getVertex(numVertices-1)-rod->getVertex(0);
	e = check->shorten;
    double e1;

    e1 = temp.norm();

	outfile<<currentTime<<" "<< f << " " <<f1<<" "<< temp1[1]<<" "<< e <<" "<<e1<<" "<<dampingF.norm()<<" "<<inertial.norm()<<endl;	
}


void world::setRodStepper()
{
	// Set up geometry
	rodGeometry();	

	// Create the rod 
	rod = new elasticRod(vertices, vertices, density, rodRadius, deltaTime,
		youngM, shearM, RodLength, theta);

	// Find out the tolerance, e.g. how small is enough?
	characteristicForce = M_PI * pow(rodRadius ,4)/4.0 * youngM / pow(RodLength, 2);
	forceTol = tol * characteristicForce;
	
	// Set up boundary condition
	rodBoundaryCondition();
	
	// setup the rod so that all the relevant variables are populated
	rod->setup();
	// End of rod setup
	
	// set up the time stepper
	stepper = new timeStepper(*rod);
	totalForce = stepper->getForce();

	// declare the forces
	m_stretchForce = new elasticStretchingForce(*rod, *stepper);
	m_bendingForce = new elasticBendingForce(*rod, *stepper);
	m_twistingForce = new elasticTwistingForce(*rod, *stepper);
	m_inertialForce = new inertialForce(*rod, *stepper);	
	m_gravityForce = new externalGravityForce(*rod, *stepper, gVector);
	m_dampingForce = new dampingForce(*rod, *stepper, viscosity, eta_per, eta_par);
    
    check = new collision(*rod, *stepper, friction);
	Nstep = totalTime/deltaTime;
    
	// Allocate every thing to prepare for the first iteration
	rod->updateTimeStep();
	
	timeStep = 0;
	currentTime = 0.0;
}

// Setup geometry
void world::rodGeometry()
{
	vertices = MatrixXd(numVertices, 3);
    

    // get the intiial configuration from the .txt file
	ifstream myfile("simDER.txt"); 

	int row1 =  numVertices; 

	MatrixXd data = MatrixXd(row1, 4);
    double a ;
	if (myfile.is_open())
	{
		for (int i = 0; i<row1* 4; i++)
		{
			myfile>>a;
			if (i%4 == 0)
				data(i/4, 0) = a;
			else if (i%4 == 1)
				data(i/4, 1) = a;
			else if (i%4 == 2)
				data(i/4, 2) = a;
			else if (i%4 == 3)
				data(i/4, 3) = a;
		}
	}
    theta = VectorXd::Zero(numVertices - 1);
	for (int i = 0; i< numVertices; i++)
	{
		vertices(i, 0) = data(i,0);
		vertices(i, 1) = data(i,1);
		vertices(i, 2) = data(i,2);

		if (i < numVertices-1)
		{
			theta(i) = data(i,3);
		}
	}
   // theta = VectorXd::Zero(numVertices - 1);

}


void world::rodBoundaryCondition()
{

	rod->setVertexBoundaryCondition(rod->getVertex(0),0);
    rod->setThetaBoundaryCondition(rod->getTheta(0),0);
    rod->setVertexBoundaryCondition(rod->getVertex(1),1);

	rod->setVertexBoundaryCondition(rod->getVertex(numVertices-1),numVertices-1);
	rod->setVertexBoundaryCondition(rod->getVertex(numVertices-2),numVertices-2);
    rod->setThetaBoundaryCondition(rod->getTheta(numVertices-2),numVertices-2);
}


void world::updateBoundary()
{
	Vector3d u;
	u(0) = 0;
	u(1)  = 0.001;
	u(2) = 0;
    
    
	if (currentTime>statictime && currentTime<= statictime+pulltime)
	{
		rod->setVertexBoundaryCondition(rod->getVertex(0)-u*deltaTime  ,0);
		rod->setVertexBoundaryCondition(rod->getVertex(1)-u*deltaTime ,1);

		rod->setVertexBoundaryCondition(rod->getVertex(numVertices-2)+u*deltaTime  ,numVertices-2);
		rod->setVertexBoundaryCondition(rod->getVertex(numVertices-1)+u*deltaTime ,numVertices-1);
        
	} 

	if (currentTime >statictime+pulltime )
	{
		timeStep = Nstep;
	}

}

void world::updateCons()
{
	rod->updateMap();
	stepper->update();
    totalForce = stepper->getForce();
}
	

void world::updateTimeStep()
{	
	bool solved = false;
	iter_c = 0;
    
	updateBoundary();
	rod->updateGuess();
    //get configuration when no collisoins    
	newtonMethod(solved); 

    //obtain the possible contact pairs during this time step
    check->setZero();
    bool flag = check->candidateTest(); 

    // compute configuration of the rod with contact
    if (flag)
    {
    	check->collisionTest();
    	while (check->checker==false)
    	{
    		solved = false;
    		rod->updateGuess();
    		check->iter = iter_c;
    		newtonMethodC(solved);
    		check->collisionTest();
    		iter_c++;
    	}
    }
    // compute the contact zone for the knots
    check->contactZone();
 
    //calculate relevant forces;
    calculateForce();

    //adapt the rod configuration;

    

	
	rod->updateTimeStep();
	cout<<rod->u.norm()<<endl;

	if (render) cout << "time: " << currentTime << " iter=" << iter << " checker ="<< check->checker<<
		"iter_c= "<<iter_c<<endl;

	// cout << "time: " << currentTime << " iter=" << iter << " checker ="<< check->checker<<
	// 	"iter_c= "<<iter_c<<endl;

	currentTime += deltaTime;
		
	timeStep++;
	
	if (solved == false)
	{
		timeStep = Nstep; // we are exiting
	}
}

void world::calculateForce()
{
	stepper->setZero();

	m_inertialForce->computeFi();
    
    inertial[0] = 0;
    inertial[1] = 0;
    inertial[2] = 0;
	for (int i =0; i< rod->nv; i++ )
	{
		inertial[0] = stepper->force[4*i] + inertial[0];
		inertial[1] = stepper->force[4*i+1] + inertial[1];
		inertial[2] = stepper->force[4*i+2] + inertial[2];
	}

	dampingF[0] = 0;
	dampingF[1] = 0;
	dampingF[2] = 0;

	stepper->setZero();

	m_dampingForce->computeFd();
	for (int i =0; i< rod->nv; i++ )
	{
		dampingF[0] = stepper->force[4*i] + dampingF[0];
		dampingF[1] = stepper->force[4*i+1] + dampingF[1];
		dampingF[2] = stepper->force[4*i+2] + dampingF[2];
	}

 
    stepper->setZero();

	m_inertialForce->computeFi();
	m_stretchForce->computeFs();		
	m_bendingForce->computeFb();
	m_twistingForce->computeFt();
	m_gravityForce->computeFg();
	m_dampingForce->computeFd();

	temp[0] =stepper->force[0]+stepper->force[4];
	temp[1] = stepper->force[1]+stepper->force[5];
	temp[2] = stepper->force[2]+stepper->force[6];


	temp1[0] = stepper->force[rod->ndof-3]+stepper->force[rod->ndof-7];
	temp1[1] = stepper->force[rod->ndof-2]+stepper->force[rod->ndof-6];
	temp1[2] = stepper->force[rod->ndof-1]+stepper->force[rod->ndof-5];

	gravity[0] = 0;
	gravity[1] = 0;
	gravity[2] = -rod->rho * M_PI* pow(rod->rodRadius, 2)*rod->rodLength * 10;
	
	// cout<<temp.norm()<<" "<<temp1.norm()<<" "<<inertial.norm()<<dampingF.norm()<<endl;
    // cout<< -temp - temp1 + inertial + dampingF<<endl;
}



void world::newtonMethodC(bool &solved)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	iter = 0;
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();
			
		m_bendingForce->computeFb();
		m_bendingForce->computeJb();
		
		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();
		
		m_dampingForce->computeFd();
		m_dampingForce->computeJd();
        
        check->computeFc();
       
		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}
		normf = sqrt(normf);

		// cout<<stepper->force<<endl;
		// cout<<"iter: "<<iter<<"normf: "<<normf<<" "<<check->checker<<endl;
		if (iter == 0) 
		{
			normf0 = normf;
		}
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}
		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			exit(0);
			break;
		}
	}
}

void world::newtonMethod(bool &solved)
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	iter = 0;
	while (solved == false)
	{
		rod->prepareForIteration();
		
		stepper->setZero();

		// Compute the forces and the jacobians
		m_inertialForce->computeFi();
		m_inertialForce->computeJi();
			
		m_stretchForce->computeFs();
		m_stretchForce->computeJs();
			
		m_bendingForce->computeFb();
		m_bendingForce->computeJb();
		
		m_twistingForce->computeFt();
		m_twistingForce->computeJt();

		m_gravityForce->computeFg();
		m_gravityForce->computeJg();
		
		m_dampingForce->computeFd();
		m_dampingForce->computeJd();

		// stepper->monitor();
       
		// Compute norm of the force equations.
		normf = 0;
		for (int i=0; i < rod->uncons; i++)
		{
			normf += totalForce[i] * totalForce[i];
		}
		normf = sqrt(normf);
		// cout<<"iter: "<<iter<<"normf: "<<normf<<" "<<check->checker<<endl;
		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol )
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			rod->updateNewtonX(totalForce); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}
}

int world::simulationRunning()
{
	if (timeStep<Nstep) 
		return 1;
	else 
	{
		return -1;
	}
}

int world::numPoints()
{
	return rod->nv;
}

double world::getScaledCoordinate(int i)
{
	return rod->x[i] / RodLength;
}

double world::getCurrentTime()
{
	return currentTime;
}

double world::getTotalTime()
{
	return totalTime;
}
