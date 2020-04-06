#include <iostream>
#include <fstream>
#include <string>
#include "eigenIncludes.h"


int main()
{
	Vector3d u, v, w;
	Vector3d p0, p1, q0, q1;
	double a, b ,c ,d, e, D;
	double sc, tc, td, tn, sn , sd;

	double delta_l = 0.001;
	double SMALL_NUMBER = 1e-8;
    p0 = Vector3d(0,0,0);
    p1 = Vector3d(1,0,0);
    q0 = Vector3d(0,0,1);
    q1 = Vector3d(1,0,1);

	u = p1 - p0;
    v = q1 - q0;
    w = p0 - q0;

    u = u/delta_l;
    v = v/delta_l;
    w = w/delta_l;

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

    Vector3d p;

    u = p1 - p0;
    v = q1 - q0;
    w = p0 - q0;

    p = w + sc * u - tc * v;

    cout<<p<<endl;
    cout<<"sc: "<<sc<<" "<<"tc: "<<tc<<endl;
    cout<<p.norm()<<endl;

	return 0;

}