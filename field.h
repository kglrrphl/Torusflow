#ifndef FIELD_H
#define FIELD_H


#include <cmath>
#include <vector>

using namespace std;

double pi_2 = 2.0*M_PI, rho = (sqrt(5)-1)/2;

//velocity field

double cxcy[2][2] = {0.000000,0.618034,-1.000000,0.000000};
vector<double> phix = {1.570796,0.000000};
vector<double> phiy = {1.570796,0.000000};

#endif 
