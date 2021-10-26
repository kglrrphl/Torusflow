#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <deque>
#include <chrono>
#include <ctime>
#include <random>

#include "functions.h"

using namespace std;
double ff;

int main(int argc, char *argv[]) 
{
    
    string path = ""; //specify path where data is saved
    if (argc > 1) {
        ff = atof(argv[1]);
    }
    else { ff = 1.5; }

    double h = 0.1, h0 = 0.05*ff; //step size for video (frame duration) and for integration 
    int frames = 400, numtr = 1000, order = 8;  //number of frames, number of tracers and order of integration
    vector<deque<double>> x(numtr), y(numtr);
    vector<deque<int>> o_x(numtr), o_y(numtr);
    
    random_device seed;
    mt19937 gen(seed());
    uniform_real_distribution<float> dist(0, 1);
    
    double r = 0.2, phi = pi_2, x0 = 0.2, y0 = 0.2, rr, pphi; //r = maximum radius of tracer cloud, (x0,y0) starting position
    
    //integrate
    for (int i = 0; i < numtr; ++i)
    {
        rr = r*dist(gen); 
        pphi = phi*dist(gen);
        
        x[i].push_back(rr*cos(pphi) + x0);
        o_x[i].push_back(0);
        y[i].push_back(rr*sin(pphi) + y0); 
        o_y[i].push_back(0);
        
        integ(x[i], y[i], o_x[i], o_y[i], h, h0, frames, order);
    }
    
    //save data
    ofstream ofs_x(path + "vid_x.txt", ofstream::out);
    ofstream ofs_y(path + "vid_y.txt", ofstream::out);
    ofstream ofs_ox(path + "vid_ox.txt", ofstream::out);
    ofstream ofs_oy(path + "vid_oy.txt", ofstream::out);
    ofstream ofs_t(path + "vid_t.txt", ofstream::out);
    for (int i = 0; i < numtr; i++)
    {
        for (int j = 0; j < frames; ++j)
        {
            if (j < frames - 1){
                ofs_x << x[i][j] << ",";
                ofs_y << y[i][j] << ",";
                ofs_ox << o_x[i][j] << ",";
                ofs_oy << o_y[i][j] << ","; }
            else{
                ofs_x << x[i][j] << ";";
                ofs_y << y[i][j] << ";";
                ofs_ox << o_x[i][j] << ";";
                ofs_oy << o_y[i][j] << ";";
                ofs_t << i*h << ";" << endl;}
        }
        
        ofs_x << endl;
        ofs_y << endl;
        ofs_ox << endl;
        ofs_oy << endl;
    }
    ofs_x.close();
    ofs_y.close();
    ofs_ox.close();
    ofs_oy.close();
    ofs_t.close();
    
    cout << "all done";
    return 0;
}

