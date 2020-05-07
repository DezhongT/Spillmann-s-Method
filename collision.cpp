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


double collision::log_func(double x, double k, double c)
{
    return 1/(1+exp(-k * (x - c)));
}

double collision::boxcar_func(double x, double k)
{
    double step1 = log_func(x, k, 0.03);
    double step2 = log_func(x, k, 0.97);

    return step1 - step2;
}

double collision::fixbound(double num)
{
	if (num < 0)
		num = 0;
	else if (num > 1)
		num = 1;

	return num;
}


void collision::computeDistance(double k1, double k2)
{
    u = p1 - p0; //d1
    v = q1 - q0; //d2
    w = q0 - p0; //d12

    a = u.dot(u); //D1
    b = u.dot(v); //R
    c = v.dot(v); //D2
    d = u.dot(w); //S1
    e = v.dot(w); //S2
    D = a*c - b*b; //den
    
    double uf;
    if (D == 0)
    {
    	sc = 0;
    	tc = -e/c;
    	uf = fixbound(tc);

    	if (uf != tc)
    	{
    		sc = (uf * b + d)/a;
    		sc = fixbound(sc);
    		tc = uf;
    	}
    }
    else
    {
    	sc = (d*c - e *b)/D;
    	sc = fixbound(sc);
    	tc = (sc * b - e)/c;
    	uf = fixbound(tc);

    	if (uf != tc)
    	{
    		sc = (uf * b + d)/a;
    		sc = fixbound(sc);
    		tc = uf;
    	}
    }

    detectionP.P = u*sc-tc*v-w;
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
    vector<int> contactI;
    vector<int> contactJ;
    contactI.clear();

    //analysis the contact force 
    Vector3d temp;
    for (int i = 0; i< rod->nv; i++)
    {
        temp  = temp_fv.segment<3>(4*i);
        // cout<<"norm: "<<temp.norm()<<endl;
        if (temp.norm()!=0)
            contactI.push_back(i);
    }
    cout<<"size: "<<contactI.size()<<endl;
    if (contactI.size()!= 0)
    {
        sort(contactI.begin(), contactI.end());
        piecewise.clear();
        int s;
        s = contactI[0];

        for (int i = 0; i < contactI.size()-1; i++)
        {
            int i1 = contactI[i];
            int j1 = contactI[i+1];
            // cout<<"i: "<<i1<<endl;

            if (j1 - i1 > 20)//threshold
            {
                piece.s = s;
                piece.e = i1;
                piecewise.push_back(piece);
                s = j1;
            }
        }

        piece.s = s;
        piece.e = contactI[contactI.size()-1];
        piecewise.push_back(piece);

        // for (int i = 0; i< piecewise.size(); i++)
        // {
        //     piece = piecewise[i];
        //     cout<<"s: "<<piece.s<<" e: "<<piece.e<<endl;
        // }

        int s_braid_1, e_braid_1;

        s_braid_1 = piecewise[0].s;
        e_braid_1 = piecewise[0].e;
        braid = 0;
        for (int i = s_braid_1; i<e_braid_1+1; i++)
        {
            temp = rod->getVertex(i+1) - rod->getVertex(i);
            braid = braid + temp.norm();
        }

        if (piecewise.size()>2)
        {
            torsion = 1;
        }
        else
        {
            torsion = 0;
        }

        int loop_s, loop_e;
        loop_s = piecewise[0].e;
        loop_e = piecewise[1].s;
        loop = 0;
        for (int i = loop_s; i< loop_e; i++)
        {
            temp = rod->getVertex(i+1) - rod->getVertex(i);
            loop = loop + temp.norm();
        }
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
            
            relativeV =  (1- detectionP.sc) *rod->getVelocity(i) + detectionP.sc *rod->getVelocity(i+1) - (1-detectionP.tc)*rod->getVelocity(j) - detectionP.tc * rod->getVelocity(j+1) ;
            //obtain the tanget velocity
            relativeV = relativeV - relativeV.dot(detectionP.P)*detectionP.P;
            
            //here we did not use the multiplier of 2-iter, I think ignore it can give us a better performace
            if (relativeV.norm() > 1e-6)
            {
                 relativeV = relativeV/relativeV.norm();
                 // temp_f.block<3,1>(0,0) = temp_f.block<3,1>(0,0) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* (1- detectionP.sc) * friction ;
                 // temp_f.block<3,1>(0,1) = temp_f.block<3,1>(0,1) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* detectionP.sc * friction  ;
                 // temp_f.block<3,1>(0,2) = temp_f.block<3,1>(0,2) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)* (1- detectionP.tc) * friction ;
                 // temp_f.block<3,1>(0,3) = temp_f.block<3,1>(0,3) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)*  detectionP.tc * friction  ;

                 temp_f.block<3,1>(0,0) = temp_f.block<3,1>(0,0) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* (1- detectionP.sc) * friction * pow(0.5, iter);
                 temp_f.block<3,1>(0,1) = temp_f.block<3,1>(0,1) + relativeV * detectionP.kesi * (d0 - detectionP.distance)* detectionP.sc * friction * pow(0.5, iter) ;
                 temp_f.block<3,1>(0,2) = temp_f.block<3,1>(0,2) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)* (1- detectionP.tc) * friction * pow(0.5, iter);
                 temp_f.block<3,1>(0,3) = temp_f.block<3,1>(0,3) - relativeV * (1 - detectionP.kesi) * (d0 - detectionP.distance)*  detectionP.tc * friction * pow(0.5, iter) ;
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





