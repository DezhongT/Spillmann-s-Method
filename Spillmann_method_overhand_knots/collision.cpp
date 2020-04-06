#include "collision.h"


collision::collision(elasticRod &m_rod, timeStepper &m_stepper, double m_friction)
{
	rod = &m_rod;
	stepper = &m_stepper;
	d0 = 2*rod->rodRadius;
	delta_l = rod->rodLength/(rod->nv - 1);
    SMALL_NUMBER = 1e-10;
    interval = 2*d0/delta_l + 1;

    temp_fv = VectorXd::Zero(rod->ndof); //add force
    temp_f = MatrixXd::Zero(4,4); //add local force
    checker = false;
    
    friction = m_friction;
    iter = 0;
}

collision::~collision()
{
	;
}


double collision::H(double t)
{
    double t1;

    double k = 100;
    t1 = 1/(1+ exp(-k*(t-0.5)));

    return t1;
}

double collision::B(double t)
{
    double t1;
    
    double k = 100;
    t1 = 1/(1+exp(-k*t)) - 1/(1+exp(-k*(t-1)));

    return t1;
}


void collision::computeDistance()
{
    u = p1 - p0;
    v = q1 - q0;
    w = p0 - q0;

    a = u.dot(u);
    b = u.dot(v);
    c = v.dot(v);
    d = u.dot(w);
    e = v.dot(w);
    D = a*c - b*b;

    sc =D;
    sd =D;
    sn =D;
    td =D;
    tn =D;
    tc =D;

    if (D < SMALL_NUMBER) { // the lines are almost parallel
        sn = 0.0;         // force using point P0 on segment S1
        sd = 1.0;         // to prevent possible division by 0.0 later
        tn = e;
        td = c;
    }
    else {                 // get the closest points on the infinite lines
        sn = (b*e - c*d);
        tn = (a*e - b*d);
        if (sn < 0.0) {        // sc < 0 => the s=0 edge is visible
            sn = 0.0;
            tn = e;
            td = c;
        }
        else if (sn > sd) {  // sc > 1  => the s=1 edge is visible
            sn = sd;
            tn = e + b;
            td = c;
        }
    }

    if (tn < 0.0) {            // tc < 0 => the t=0 edge is visible
        tn = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sn = 0.0;
        else if (-d > a)
            sn = sd;
        else {
            sn = -d;
            sd = a;
        }
    }
    else if (tn > td) {      // tc > 1  => the t=1 edge is visible
        tn = td;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sn = 0;
        else if ((-d + b) > a)
            sn = sd;
        else {
            sn = (-d +  b);
            sd = a;
        }
    }

    sc = (abs(sn) < SMALL_NUMBER ? 0 : sn/ sd);
    tc = (abs(tn) < SMALL_NUMBER ? 0 : tn / td);
    
    detectionP.P =  w + sc * u - tc * v; // xi - xj

    detectionP.distance = detectionP.P.norm();

    detectionP.P = detectionP.P/detectionP.distance;
    detectionP.sc = sc;
    detectionP.tc = tc;
}

void collision::setZero()
{
	temp_fv = VectorXd::Zero(rod->ndof);
    possibleC.clear();
}

// get the possible contact paris, calcuate the contact force from those 
// possible contact pairs during each time step
bool collision::candidateTest() 
{
    // initializition
    checker = false;
    // constraint candidates selection
    double min = 1;
    double sc_1;
    double tc_1;
    for (int i = 0; i < rod->nv-3; i++)
    {
        p0 = rod->getVertex(i);
        p1 = rod->getVertex(i+1);
        for (int j = i+2; j < rod->nv-1; j++)
        {
            q0 = rod->getVertex(j);
            q1 = rod->getVertex(j+1);
            computeDistance();
            if (detectionP.distance < min)
            {
                min = detectionP.distance;
                sc_1 = detectionP.sc;
                tc_1 = detectionP.tc;
            }

            if (detectionP.distance < 1.5*d0  &&  abs(i - j) > interval)
            {
                c_p.i = i;
                c_p.j = j;
                possibleC.push_back(c_p);
            }
        }
    }
   
    cout<<"min: "<<min<<" "<<"sc: "<<sc_1<<" tc: "<<tc_1<<endl;
    if (possibleC.size()!=0)
    {
        mc_f = MatrixXd::Zero(possibleC.size(), 5); //d0 ,i , j, lamda, f_lamda
        vector<candidateP>::iterator it;
        int c = 0;
        for (it = possibleC.begin(); it != possibleC.end(); ++ it)
        {
            mc_f(c, 1) = it->i;
            mc_f(c, 2) = it->j;
            c = c+ 1;
        }
        return true;
    }
    else
        return false;
}

// get the total length of the loop part
void collision::contactZone()
{
    int c_up = 0, c_down = 0;
    if (possibleC.size()>0)
    {
    for (int i =0; i < possibleC.size(); i++)
    {
        if (mc_f(i,3)!=0)
        {
            c_up = i;
            break;
        }
    }
    int i1 = 0;
    for (int i =0; i < possibleC.size(); i++)
    {
        i1 = possibleC.size() - i - 1;
        if (mc_f(i1,3)!=0)
        {
            c_down = i1;
            break;
        }
    }
    int a1, a2 , b1 ,b2;
    a1 = mc_f(c_up, 1);
    a2 = mc_f(c_down, 1);
    b1 = mc_f(c_up, 2);
    b2 = mc_f(c_down, 2);

    shorten = 0;

    Vector3d temp;
    for (int i = a1; i< b1; i++)
    {
        temp = rod->getVertex(i) - rod->getVertex(i+1);
        shorten = shorten + temp.norm();
    }
    }
    else
    {
        shorten = 0;
    }
}

// get the contact force for each iteration during the timestep
void collision::collisionTest()
{
	checker = true;

    vector<candidateP>::iterator it;
    int c = 0;
    for (it = possibleC.begin(); it != possibleC.end(); ++ it)
    {
        int i = it->i;
        int j = it->j;
        p0 = rod->getVertex(i);
        p1 = rod->getVertex(i+1);
        q0 = rod->getVertex(j);
        q1 = rod->getVertex(j+1);
        computeDistance();
        mc_f(c, 0) = detectionP.distance - d0;
        if ( d0 - detectionP.distance > 1e-4 * d0)
        {
        	mc_f(c, 3) = 1;
        	checker = false;
        	detectionP.kesi = (rod->massArray(4*j)*(1 - detectionP.tc) + rod->massArray(4*(j+1))
                              * detectionP.tc)/(rod->massArray(4*i)*(1 - detectionP.sc) +
                               rod->massArray(4*(i+1))*detectionP.sc+ rod->massArray(4*j)*
                               (1-detectionP.tc) + rod->massArray(4*(j+1))*detectionP.tc); 
            
            temp_f.block<3,1>(0,0) =   detectionP.P * detectionP.kesi * (detectionP.distance - d0) * (1-detectionP.sc);
            temp_f.block<3,1>(0,1) =   detectionP.P * detectionP.kesi * (detectionP.distance -  d0) * detectionP.sc;
            temp_f.block<3,1>(0,2) =   detectionP.P* (1 - detectionP.kesi) * (d0 - detectionP.distance) * (1-detectionP.tc);
            temp_f.block<3,1>(0,3) =   detectionP.P* (1 - detectionP.kesi) * (d0 - detectionP.distance) * detectionP.tc;
            
            //add friction
            relativeV =  (1- detectionP.sc) *rod->getVelocity(i) + detectionP.sc *rod->getVelocity(i+1) - (1-detectionP.tc)*rod->getVelocity(j) - detectionP.tc * rod->getVelocity(j+1) ;
            //obtain the tanget velocity
            relativeV = relativeV - relativeV.dot(detectionP.P)*detectionP.P;
            
            //here we did not use the multiplier of 2-iter, I think ignore it can give us a better performace
            if (relativeV.norm() > 1e-6)
            {
                 relativeV = relativeV/relativeV.norm();
                 temp_f.block<3,1>(0,0) = temp_f.block<3,1>(0,0) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* (1- detectionP.sc) * friction ;
                 temp_f.block<3,1>(0,1) = temp_f.block<3,1>(0,1) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* detectionP.sc * friction  ;
                 temp_f.block<3,1>(0,2) = temp_f.block<3,1>(0,2) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)* (1- detectionP.tc) * friction ;
                 temp_f.block<3,1>(0,3) = temp_f.block<3,1>(0,3) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)*  detectionP.tc * friction  ;

                 // temp_f.block<3,1>(0,0) = temp_f.block<3,1>(0,0) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* (1- detectionP.sc) * friction * pow(0.5, iter);
                 // temp_f.block<3,1>(0,1) = temp_f.block<3,1>(0,1) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* detectionP.sc * friction * pow(0.5, iter) ;
                 // temp_f.block<3,1>(0,2) = temp_f.block<3,1>(0,2) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)* (1- detectionP.tc) * friction * pow(0.5, iter);
                 // temp_f.block<3,1>(0,3) = temp_f.block<3,1>(0,3) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)*  detectionP.tc * friction * pow(0.5, iter) ;
            }
            temp_f(3, 0) = i;
            temp_f(3, 1) = i + 1;
            temp_f(3, 2) = j;
            temp_f(3, 3) = j + 1;

            for (int c1 = 0; c1 < 4; c1 ++)
            {
                temp_f.block<3,1>(0,c1) = temp_f.block<3,1>(0,c1) *rod->massArray(4*i) / pow(rod->dt, 2);
                for (int c2 = 0; c2 < 3; c2 ++)
                {
                     temp_fv(4*temp_f(3,c1)+c2) = temp_fv(4*temp_f(3,c1)+c2) + temp_f(c2, c1);
                }
            }
        }
        c = c + 1;
    }
}


void collision::monitor()
{
	cout<<temp_fv<<endl;
}



void collision::computeFc()
{
	for (int i =0 ; i<rod->ndof; i++)
	{
		stepper->addForce(i, temp_fv(i));
	}
}





