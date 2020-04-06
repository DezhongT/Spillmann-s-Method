#include <iostream>
#include <fstream>
#include <string>
#include "eigenIncludes.h"
// #include <MATH>


double H(double t)
{
    double t1;

    double k = 100;
    t1 = 1/(1+ exp(-k*(t-0.5)));

    return t1;
}

double B(double t)
{
    double t1;
    
    double k = 100;
    t1 = 1/(1+exp(-k*t)) - 1/(1+exp(-k*(t-1)));

    return t1;
}


int main()
{
	Vector3d u, v, w;
	Vector3d p0, p1, q0, q1;
    Vector3d p;
	double a, b ,c ,d, e, D;
	double sc, tc, td, tn, sn , sd;

	double delta_l = 0.001;
	double SMALL_NUMBER = 1e-8;
    p0 = Vector3d(0,0,0);//p1s
    p1 = Vector3d(1,0,0); //p1e
    q0 = Vector3d(-0.2,0,1); //q1s
    q1 = Vector3d(0.8,0,1); //q1e

	u = p1 - p0; //d1
    v = q1 - q0; //d2
    w = q0 - p0; //d12

    a = u.dot(u); //D1
    b = u.dot(v); //R
    c = v.dot(v); //D2
    d = u.dot(w); //S1
    e = v.dot(w); //S2
    D = a*c - b*b; //den
    

    D = D + 1e-9;
    // if (D == 0)
    //     D = 1e-8;

    //t is sc, u is tc;

    sc = (c * d - b*e)/D;

    sc = H(sc);

    tc = (sc * b - e)/c;

    sc = (1 - B(tc)) * (tc * b + d)/a + B(tc) * sc;

    tc = H(tc);

    // cout<<"tc: "<<tc<" sc: "<<sc<<endl;

    p =  sc * u - tc * v - w;

    cout<<p<<endl;
    cout<<"sc: "<<sc<<" "<<"tc: "<<tc<<endl;
    cout<<p.norm()<<endl;

	return 0;

}