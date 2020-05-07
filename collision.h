#ifndef COLLISION_H
#define COLLISION_H

#include "eigenIncludes.h"
#include "elasticRod.h"
#include "timeStepper.h"


struct collisionP{
    int node1;
    int node2;
    Vector3d P; //P is the distance vector between two segments
    double distance;

    double sc; //ratio in segment i -> i+1
    double tc; //ratio in segment j -> j+1
    double kesi;
};

struct candidateP{
    int i;
    int j;
};

struct node
{
    int s;
    int e;
};

class collision
{
public:
	collision(elasticRod &m_rod, timeStepper &m_stepper, double m_friction);
	~collision();
	elasticRod *rod;
	timeStepper *stepper;

	void collisionTest();
	// void clear();

	bool checkcollision();

    void update();

    void prepareforcheck();

    bool checker;
    //check whether collision happens
    bool hold;

    void setZero();

    bool candidateTest();
    vector<int> holdforTying();

    VectorXd temp_fv;
    MatrixXd mc0;
    void contactZone();

    MatrixXd mc_f;

    double shorten;

    int iter;
    void computeFc();

    void monitor();

    double fixbound(double num);

    double braid;
    double loop;

    double velocity;
    double torsion;

private:
	// parameter for collision check
	double SMALL_NUMBER;
	Vector3d p0;
	Vector3d p1;
	Vector3d q0;
	Vector3d q1;
    Vector3d u;
    Vector3d v;
    Vector3d w;
    double tn;
    double td;
    double sn;
    double sd;
    double tc;
    double sc;
    double a;
    double b;
    double c;
    double d;
    double e;
    double D;
    collisionP detectionP;
    int interval;


    //mininum distance between rod segment
    double d0;
    //ref length between two nodes;
    double delta_l;


    //Matrix for distance
    MatrixXd md;

    MatrixXd mv;

    MatrixXd mv0;

    MatrixXd mc;

    //vector for collision candidate
    candidateP c_p;
    vector <candidateP> possibleC;

    void computeDistance(double k1 = 6.5, double k2 = 10);

    double penetrationL;

    Vector3d temp_v;
    Vector3d temp_v0;

    MatrixXd temp_f;

    double friction;
    // double v_constant;

    Vector3d relativeV;

    double log_func(double x, double k, double c = 0.5);
    double boxcar_func(double x, double k = 100);

    
    vector<node> piecewise;
    node piece;


};
#endif