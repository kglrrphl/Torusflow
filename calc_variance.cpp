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


#include "functions.h"

using namespace std;
double ff;

int main(int argc, char *argv[])
{   
    clock_t t, tt;
    string path = ""; //specify path where data is saved
    if (argc > 1) {
        ff = atof(argv[1]);
    }
    else { ff = 1.5; }
    
    //simulation parameters
    double xstart = 1.5, ystart = 1.5; //starting values
    double h = 1.0*ff, h0 = 0.05*ff; //sampling interval of trajectory and initial step size
    int endtime_zeta = 1e6, runs = 1, order = 6, ntp = 100; // endtime of zeta, number of runs, order of integration, number of sampling points
    
    //declaration of variables
    deque<double> x(1, xstart), y(1, ystart), vy; //trajectories on torus
    deque<int> o_x(1), o_y(1); //revolutions
    double psi_0, psi_end;

    double v_mean = 0, v_mean_o = 0, T_mean = 0, v_mean_end, p_nc;
    double mxi, sxi, fact; //auxiliary variables for averaging
    vector<double> xi_mean(ntp), xi_var(ntp); //mean and variance at sample points
    vector<double> drift(ntp); //estimated drift at sample points
    vector<int> tp(ntp); //sampling points of zeta

    //infos for dong jobs
    time_t datum = chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << ctime(&datum) << endl;
    cout << "amplitude = " << ff << endl << "simulation time = " << endtime_zeta << " * " << runs << endl;
    cout << "path: " << path << endl;

    //logspaced time points
    logspace(tp, 0, log10((double)endtime_zeta));
    
    write2file(tp, path + "time.txt");
    write2file(h, path + "tscale.txt");
    //write parameters to file
    ofstream ofs(path + "parameters.txt", ofstream::out);
    ofs << "f = " << ff << endl;
    ofs << "starting point: (" << xstart << ", " << ystart << ")" << endl;
    ofs << "time scale: " << h << endl;
    ofs << "endtime of zeta: " << endtime_zeta << endl;
    ofs << "simulation runs: " << runs << endl;
    ofs << "stepsize of integration: " << h0/ff << endl;
    ofs << "order of integration: " << order << endl;
    ofs.close();
    
    //streamfunction at start
    psi_0 = streamf(xstart, ystart);
    cout << "psi(0) = " << psi_0 << endl;
    
    //integrate trajectory once up front
    integ(x, y, o_x, o_y, h, h0, endtime_zeta, order, vy);
    
    //estimate mean velocity and drift
    int delay = 100;
	for (int i = 0; i < endtime_zeta - delay; i++)
	{
		v_mean += y[i + delay] - y[i];
		v_mean_o += o_y[i + delay] - o_y[i];
	}
    
    v_mean = (v_mean + pi_2*v_mean_o)/((endtime_zeta - delay)*(delay*h));
	        
    for (int i = 0; i < ntp; ++i)
    {
            drift[i] = v_mean*(tp[i]*h);
    }
    cout << "v_mean = " << v_mean << endl;
    

    t = clock();
    //loop over runs
    for (int n = 0; n < runs; ++n) 
    {
            cout << "run " << n+1 << " / " << runs << endl;
            tt = clock();

            //integration
            if (n > 0) { integ(x, y, o_x, o_y, h, h0, endtime_zeta, order); }
            else { 
                integ(x, y, o_x, o_y, h, h0, endtime_zeta, order, vy); 
                
                //mean passage time
                int ind = 0, norm = 0;
                for (int i = 1; i < o_y.size(); ++i) { 
                    if (o_y[i] > o_y[ind]) {
                        T_mean += i - ind;
                        ind = i;
                        norm += 1; }}
                T_mean /= ff*norm;
                cout << "T_mean = " << T_mean << endl;
                
                //correlation function and power spectrum
                for (int i = 0; i < vy.size(); ++i) { vy[i] -= v_mean; }   
                int zpads = (1 << int(ceil(log2(vy.size())))) - vy.size();
                vy.insert(vy.end(), zpads, 0);
                vy.insert(vy.end(), vy.size(), 0);
                
                vector<double> ps(vy.size());
                vector<double> cvv(vy.size());
                
                powerspec(vy, ps);
                correl(vy, vy, cvv);
                for (int j = 1; j < cvv.size(); ++j) { cvv[j] /= cvv[0]; }
                cvv.erase(cvv.begin() + (cvv.size() >> 1), cvv.end());
                 
                write2file(T_mean, path + "t_mean.txt");
                write2file(cvv, path + "correlation.txt");
                write2file(ps, path + "spectrum.txt"); 
            }
            
            //check for error
            psi_end = streamf(x.back() + pi_2 * o_x.back(), y.back() + pi_2 * o_y.back());
            cout << "psi(end) = " << psi_end << endl;
            cout << "error = " << abs(psi_0 - psi_end) << endl;
            
            fact = endtime_zeta*n/(n + 1);
            
            //loop over sampling points of zeta
            for (int i = 0; i < ntp; ++i) 
            {
                mxi = 0;
                sxi = 0;

                var_doublesum(mxi, sxi, y, o_y, endtime_zeta, tp[i], drift[i]);

                //update mean and variance
                if (n > 0){
                    xi_var[i] += sxi + (xi_mean[i] - mxi)*(xi_mean[i] - mxi)*fact;
                }
                else{
                    xi_var[i] = sxi;
                }

                xi_mean[i] += (mxi - xi_mean[i]) / (n + 1);
            }

            //erase obsolete values
            x.erase(x.begin(), x.begin() + endtime_zeta);
            y.erase(y.begin(), y.begin() + endtime_zeta);
            o_x.erase(o_x.begin(), o_x.begin() + endtime_zeta);
            o_y.erase(o_y.begin(), o_y.begin() + endtime_zeta);

            //reset revolutions
            for (int i = o_x.size() - 1; i >= 0; --i)
            {
                    o_x[i] -= o_x[0];
                    o_y[i] -= o_y[0];
            }
            
            psi_0 = streamf(x.back() + pi_2 * o_x.back(), y.back() + pi_2 * o_y.back());
            tt = clock() - tt;
            cout << "took " << ((float)tt) / CLOCKS_PER_SEC << "s" << endl;
    }

    for (int i = 0; i < ntp; ++i)
    {
        xi_var[i] /= runs*endtime_zeta;   
        xi_mean[i] += drift[i];
    }
    
    v_mean_end = abs(xi_mean.back())/(tp.back()*h), p_nc = 1/v_mean_end;
    
    //write results to file
    write2file(xi_var, path + "variance.txt");
    write2file(xi_mean, path + "mean.txt");
    write2file(p_nc, path + "p_nc.txt");
    
    t = clock() - t;
    t = ((float)t) / CLOCKS_PER_SEC;
    write2file(t, path + "runtime.txt");
    cout << "total runtime: " << t << "s" << endl;
    cout << "all done" << endl;

    return 0;
}
