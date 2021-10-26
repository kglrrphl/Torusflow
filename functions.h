#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <deque>

#include "field.h"

using namespace std;
extern double ff;

//stream function
double streamf(double x, double y)
{
    int ind = phix.size();
    double psi = 0;
    for (int i = 0; i < ind; ++i){
        for (int j = 0; j < ind; ++j)
        {
           psi += cxcy[i][j]*sin(i*x + phix[i])*sin(j*y + phiy[j]);
        }}
        
    psi = ff*psi + (-x + rho*y);
    
    return(psi);
}     

//integration of arbitrary order
void integ(deque<double>& x, deque<double>& y, deque<int>& o_x, deque<int>& o_y, double h, double h0, int steps, int order)
{
    double xx = x.back(), yy = y.back(), t_elapsed, hh = h0; 
    int oo_x = o_x.back(), oo_y = o_y.back(), count = 1, harm = phix.size(), skiprow[harm], skipcol[harm];
    double a = rho, b = -1.0, kovern[order][harm], hpownovern[order], cxxcy[harm][harm], cxcyy[harm][harm];
    double dx[order], dy[order], dx1[order][harm], dx2[order][harm], dy1[order][harm], dy2[order][harm];
    
    
    for (int n = 0; n < order; ++n){
        for (int k = 0; k < harm; ++k)
        {
            if (n == 0){ kovern[n][k] = 1; }
            else{ kovern[n][k] = double(k)/double(n); }
        }}
        
    for (int i = 0; i < harm; ++i){
        for (int j = 0; j < harm; ++j)
        { 
            cxcyy[i][j] = cxcy[i][j]*j*ff; 
            cxxcy[i][j] = cxcy[i][j]*i*ff; 
            
            if (cxcy[i][j] != 0){ 
                skiprow[i] = 1;
                skipcol[j] = 1;}
        }}
    
    
    //starting values
    dx[0] = xx;
    dy[0] = yy;
    for (int k = 0; k < harm; k++)
    {
        if (skiprow[k] == 1){
            dx1[0][k] = sin(k*xx + phix[k]);
            dx2[0][k] = cos(k*xx + phix[k]);}
        if (skipcol[k] == 1){
            dy1[0][k] = sin(k*yy + phiy[k]);
            dy2[0][k] = cos(k*yy + phiy[k]);}
    }
        
    //big loop over external timesteps
    while (count <= steps)
    {
        t_elapsed = 0.0;
              
        while (t_elapsed < h)
        {
            if ((t_elapsed + hh) > h)
            {hh = h - t_elapsed;}
            
            //reset derivatives
            for (int n = 1; n < order; ++n){
                if (n == 1){
                    dx[n] = a;
                    dy[n] = -b;}
                else{
                    dx[n] = 0;
                    dy[n] = 0;}
                
            }
                
            for (int k = 0; k < harm; ++k){
                if (skiprow[k] == 1){
                    for (int n = 1; n < order; ++n){
                        dx1[n][k] = 0;
                        dx2[n][k] = 0;}}
                if (skipcol[k] == 1){
                    for (int n = 1; n < order; ++n){
                        dy1[n][k] = 0;
                        dy2[n][k] = 0;}}
            }
            
            //calculate derivatives
            for (int n = 1; n < order; ++n){
                for (int m = 0; m < n; ++m){
                    for (int j = 0; j < harm; ++j){
                        for (int i = 0; i < harm; ++i)
                        {  
                            if (cxcy[i][j] == 0){ continue; }
                            dx[n] += cxcyy[i][j]*dx1[n - 1 - m][i]*dy2[m][j];
                            dy[n] -= cxxcy[i][j]*dy1[n - 1 - m][j]*dx2[m][i];
                        }}}
                  
                for (int k = 1; k < harm; ++k){ 
                    if (skiprow[k] == 1){
                        for (int m = 0; m < n; ++m){   
                            dx1[n][k] += dx2[n - 1 - m][k]*dx[m + 1];
                            dx2[n][k] += dx1[n - 1 - m][k]*dx[m + 1];}}
                    if (skipcol[k] == 1){
                        for (int m = 0; m < n; ++m){
                            dy1[n][k] += dy2[n - 1 - m][k]*dy[m + 1];
                            dy2[n][k] += dy1[n - 1 - m][k]*dy[m + 1];}}
                    
                    dx1[n][k] *= kovern[n][k];
                    dx2[n][k] *= -kovern[n][k];
                    dy1[n][k] *= kovern[n][k];
                    dy2[n][k] *= -kovern[n][k];
                }
            }
            
            //new values
            for (int n = 1; n < order; ++n){hpownovern[n] = pow(hh, n)/n;}
            
            for (int n = 1; n < order; ++n)
            {
                dx[0] += hpownovern[n]*dx[n];
                dy[0] += hpownovern[n]*dy[n];
            }
            
            t_elapsed += hh;
            
            //map onto torus
            if (dx[0] >= pi_2){
                dx[0] -= pi_2;
                oo_x += 1;}
            else if (dx[0] < 0){
                dx[0] += pi_2;
                oo_x -= 1;
            }
            
            if (dy[0] >= pi_2){
                dy[0] -= pi_2;
                oo_y += 1;}
            else if (dy[0] < 0){
                dy[0] += pi_2;
                oo_y -= 1;
            } 
            
            xx = dx[0];
            yy = dy[0];
            
            for (int k = 1; k < harm; k++)
            {
                if (skiprow[k] == 1){
                    dx1[0][k] = sin(k*xx + phix[k]);
                    dx2[0][k] = cos(k*xx + phix[k]);}
                if (skipcol[k] == 1){
                    dy1[0][k] = sin(k*yy + phiy[k]);
                    dy2[0][k] = cos(k*yy + phiy[k]);}
            }
        }
        
        x.push_back(xx);
        o_x.push_back(oo_x);
        y.push_back(yy);
        o_y.push_back(oo_y);
        
        count += 1;
        
        hh = h0;
    }
}

//integration overload for velocity in y direction
void integ(deque<double>& x, deque<double>& y, deque<int>& o_x, deque<int>& o_y, double h, double h0, int steps, int order, deque<double>& vy)
{
    double xx = x.back(), yy = y.back(), t_elapsed, hh = h0; 
    int oo_x = o_x.back(), oo_y = o_y.back(), count = 1, harm = phix.size(), skiprow[harm], skipcol[harm];
    double a = rho, b = -1.0, kovern[order][harm], hpownovern[order], cxxcy[harm][harm], cxcyy[harm][harm];
    double dx[order], dy[order], dx1[order][harm], dx2[order][harm], dy1[order][harm], dy2[order][harm];
    
    
    for (int n = 0; n < order; ++n){
        for (int k = 0; k < harm; ++k)
        {
            if (n == 0){ kovern[n][k] = 1; }
            else{ kovern[n][k] = double(k)/double(n); }
        }}
        
    for (int i = 0; i < harm; ++i){
        for (int j = 0; j < harm; ++j)
        { 
            cxcyy[i][j] = cxcy[i][j]*j*ff; 
            cxxcy[i][j] = cxcy[i][j]*i*ff; 
            
            if (cxcy[i][j] != 0){ 
                skiprow[i] = 1;
                skipcol[j] = 1;}
        }}
    
    
    //starting values
    dx[0] = xx;
    dy[0] = yy;
    for (int k = 0; k < harm; k++)
    {
        if (skiprow[k] == 1){
            dx1[0][k] = sin(k*xx + phix[k]);
            dx2[0][k] = cos(k*xx + phix[k]);}
        if (skipcol[k] == 1){
            dy1[0][k] = sin(k*yy + phiy[k]);
            dy2[0][k] = cos(k*yy + phiy[k]);}
    }
        
    //big loop over external timesteps
    while (count <= steps)
    {
        t_elapsed = 0.0;
              
        while (t_elapsed < h)
        {
            if ((t_elapsed + hh) > h)
            {hh = h - t_elapsed;}
            
            //reset derivatives
            for (int n = 1; n < order; ++n){
                if (n == 1){
                    dx[n] = a;
                    dy[n] = -b;}
                else{
                    dx[n] = 0;
                    dy[n] = 0;}
                
            }
                
            for (int k = 0; k < harm; ++k){
                if (skiprow[k] == 1){
                    for (int n = 1; n < order; ++n){
                        dx1[n][k] = 0;
                        dx2[n][k] = 0;}}
                if (skipcol[k] == 1){
                    for (int n = 1; n < order; ++n){
                        dy1[n][k] = 0;
                        dy2[n][k] = 0;}}
            }
            
            //calculate derivatives
            for (int n = 1; n < order; ++n){
                for (int m = 0; m < n; ++m){
                    for (int j = 0; j < harm; ++j){
                        for (int i = 0; i < harm; ++i)
                        {  
                            if (cxcy[i][j] == 0){ continue; }
                            dx[n] += cxcyy[i][j]*dx1[n - 1 - m][i]*dy2[m][j];
                            dy[n] -= cxxcy[i][j]*dy1[n - 1 - m][j]*dx2[m][i];
                        }}}
                  
                for (int k = 1; k < harm; ++k){ 
                    if (skiprow[k] == 1){
                        for (int m = 0; m < n; ++m){   
                            dx1[n][k] += dx2[n - 1 - m][k]*dx[m + 1];
                            dx2[n][k] += dx1[n - 1 - m][k]*dx[m + 1];}}
                    if (skipcol[k] == 1){
                        for (int m = 0; m < n; ++m){
                            dy1[n][k] += dy2[n - 1 - m][k]*dy[m + 1];
                            dy2[n][k] += dy1[n - 1 - m][k]*dy[m + 1];}}
                    
                    dx1[n][k] *= kovern[n][k];
                    dx2[n][k] *= -kovern[n][k];
                    dy1[n][k] *= kovern[n][k];
                    dy2[n][k] *= -kovern[n][k];
                }
            }
            
            //new values
            for (int n = 1; n < order; ++n){hpownovern[n] = pow(hh, n)/n;}
            
            for (int n = 1; n < order; ++n)
            {
                dx[0] += hpownovern[n]*dx[n];
                dy[0] += hpownovern[n]*dy[n];
            }
            
            t_elapsed += hh;
            
            //map onto torus
            if (dx[0] >= pi_2){
                dx[0] -= pi_2;
                oo_x += 1;}
            else if (dx[0] < 0){
                dx[0] += pi_2;
                oo_x -= 1;
            }
            
            if (dy[0] >= pi_2){
                dy[0] -= pi_2;
                oo_y += 1;}
            else if (dy[0] < 0){
                dy[0] += pi_2;
                oo_y -= 1;
            } 
            
            xx = dx[0];
            yy = dy[0];
            
            for (int k = 1; k < harm; k++)
            {
                if (skiprow[k] == 1){
                    dx1[0][k] = sin(k*xx + phix[k]);
                    dx2[0][k] = cos(k*xx + phix[k]);}
                if (skipcol[k] == 1){
                    dy1[0][k] = sin(k*yy + phiy[k]);
                    dy2[0][k] = cos(k*yy + phiy[k]);}
            }
        }
        
        x.push_back(xx);
        o_x.push_back(oo_x);
        y.push_back(yy);
        o_y.push_back(oo_y);
        vy.push_back(dy[1]);
        
        count += 1;
        
        hh = h0;
    }
}

//integration overload for velocity in x & y direction
void integ(deque<double>& x, deque<double>& y, deque<int>& o_x, deque<int>& o_y, double h, double h0, int steps, int order, deque<double>& vx, deque<double>& vy)
{
    double xx = x.back(), yy = y.back(), t_elapsed, hh = h0; 
    int oo_x = o_x.back(), oo_y = o_y.back(), count = 1, harm = phix.size(), skiprow[harm], skipcol[harm];
    double a = rho, b = -1.0, kovern[order][harm], hpownovern[order], cxxcy[harm][harm], cxcyy[harm][harm];
    double dx[order], dy[order], dx1[order][harm], dx2[order][harm], dy1[order][harm], dy2[order][harm];
    
    
    for (int n = 0; n < order; ++n){
        for (int k = 0; k < harm; ++k)
        {
            if (n == 0){ kovern[n][k] = 1; }
            else{ kovern[n][k] = double(k)/double(n); }
        }}
        
    for (int i = 0; i < harm; ++i){
        for (int j = 0; j < harm; ++j)
        { 
            cxcyy[i][j] = cxcy[i][j]*j*ff; 
            cxxcy[i][j] = cxcy[i][j]*i*ff; 
            
            if (cxcy[i][j] != 0){ 
                skiprow[i] = 1;
                skipcol[j] = 1;}
        }}
    
    
    //starting values
    dx[0] = xx;
    dy[0] = yy;
    for (int k = 0; k < harm; k++)
    {
        if (skiprow[k] == 1){
            dx1[0][k] = sin(k*xx + phix[k]);
            dx2[0][k] = cos(k*xx + phix[k]);}
        if (skipcol[k] == 1){
            dy1[0][k] = sin(k*yy + phiy[k]);
            dy2[0][k] = cos(k*yy + phiy[k]);}
    }
        
    //big loop over external timesteps
    while (count <= steps)
    {
        t_elapsed = 0.0;
              
        while (t_elapsed < h)
        {
            if ((t_elapsed + hh) > h)
            {hh = h - t_elapsed;}
            
            //reset derivatives
            for (int n = 1; n < order; ++n){
                if (n == 1){
                    dx[n] = a;
                    dy[n] = -b;}
                else{
                    dx[n] = 0;
                    dy[n] = 0;}
                
            }
                
            for (int k = 0; k < harm; ++k){
                if (skiprow[k] == 1){
                    for (int n = 1; n < order; ++n){
                        dx1[n][k] = 0;
                        dx2[n][k] = 0;}}
                if (skipcol[k] == 1){
                    for (int n = 1; n < order; ++n){
                        dy1[n][k] = 0;
                        dy2[n][k] = 0;}}
            }
            
            //calculate derivatives
            for (int n = 1; n < order; ++n){
                for (int m = 0; m < n; ++m){
                    for (int j = 0; j < harm; ++j){
                        for (int i = 0; i < harm; ++i)
                        {  
                            if (cxcy[i][j] == 0){ continue; }
                            dx[n] += cxcyy[i][j]*dx1[n - 1 - m][i]*dy2[m][j];
                            dy[n] -= cxxcy[i][j]*dy1[n - 1 - m][j]*dx2[m][i];
                        }}}
                  
                for (int k = 1; k < harm; ++k){ 
                    if (skiprow[k] == 1){
                        for (int m = 0; m < n; ++m){   
                            dx1[n][k] += dx2[n - 1 - m][k]*dx[m + 1];
                            dx2[n][k] += dx1[n - 1 - m][k]*dx[m + 1];}}
                    if (skipcol[k] == 1){
                        for (int m = 0; m < n; ++m){
                            dy1[n][k] += dy2[n - 1 - m][k]*dy[m + 1];
                            dy2[n][k] += dy1[n - 1 - m][k]*dy[m + 1];}}
                    
                    dx1[n][k] *= kovern[n][k];
                    dx2[n][k] *= -kovern[n][k];
                    dy1[n][k] *= kovern[n][k];
                    dy2[n][k] *= -kovern[n][k];
                }
            }
            
            //new values
            for (int n = 1; n < order; ++n){hpownovern[n] = pow(hh, n)/n;}
            
            for (int n = 1; n < order; ++n)
            {
                dx[0] += hpownovern[n]*dx[n];
                dy[0] += hpownovern[n]*dy[n];
            }
            
            t_elapsed += hh;
            
            //map onto torus
            if (dx[0] >= pi_2){
                dx[0] -= pi_2;
                oo_x += 1;}
            else if (dx[0] < 0){
                dx[0] += pi_2;
                oo_x -= 1;
            }
            
            if (dy[0] >= pi_2){
                dy[0] -= pi_2;
                oo_y += 1;}
            else if (dy[0] < 0){
                dy[0] += pi_2;
                oo_y -= 1;
            } 
            
            xx = dx[0];
            yy = dy[0];
            
            for (int k = 1; k < harm; k++)
            {
                if (skiprow[k] == 1){
                    dx1[0][k] = sin(k*xx + phix[k]);
                    dx2[0][k] = cos(k*xx + phix[k]);}
                if (skipcol[k] == 1){
                    dy1[0][k] = sin(k*yy + phiy[k]);
                    dy2[0][k] = cos(k*yy + phiy[k]);}
            }
        }
        
        x.push_back(xx);
        o_x.push_back(oo_x);
        y.push_back(yy);
        o_y.push_back(oo_y);
        vx.push_back(dx[1]);
        vy.push_back(dy[1]);
        
        count += 1;
        
        hh = h0;
    }
}

//logspaced vector
template <typename T>
void logspace(vector<T>& v, double a, double b)
{
	double ival = (b - a) / (((double)v.size()) - 1.0);

	for (size_t i = 0; i < v.size(); i++)
	{
            v[i] = T(pow(10, a + i * ival));
            if (i > 0 && is_same<T, int>::value){
                if (v[i] <= v[i - 1])
                {v[i] = v[i -1] + 1;}
            }
	}
}

//linspaced vector
template <typename T>
void linspace(vector<T>& v, double a, double b)
{
	double ival = (b - a) / ((double(v.size())) - 1.0);

	for (size_t i = 0; i < v.size(); i++)
	{
            v[i] = T((a + i*ival));
	}
}

//averaging routine
void var_doublesum(double& mean, double& psum, deque<double> v, deque<int> v_o, int N, int T, double d)
{
    double xi;
    
    for (int j = 0; j < N; ++j) // loop over trajectory
    {
        xi = (v[j + T] - v[j]) + pi_2*(v_o[j + T] - v_o[j]) - d;
        mean += xi;
        psum += xi*xi;
    }
    
    
    psum = psum - mean*mean/N;
    mean = mean/N;
}

//averaging routine for crossterm
void cross_doublesum(double& psum, deque<double> vx, deque<int> vx_o, deque<double> vy, deque<int> vy_o, int N, int T, double dx, double dy)
{
    double xi_x, xi_y;
    
    for (int j = 0; j < N; ++j) // loop over trajectory
    {
        xi_x = (vx[j + T] - vx[j]) + pi_2*(vx_o[j + T] - vx_o[j]) - dx;
        xi_x = (vy[j + T] - vy[j]) + pi_2*(vy_o[j + T] - vy_o[j]) - dy;
        psum += xi_x*xi_y;
    }
}

//write deque to file
template <typename T>
void write2file(deque<T> d, string filename)
{
	ofstream ofs(filename, ofstream::out);

	for (size_t i = 0; i < d.size(); i++)
	{
		ofs << d[i] << ";" << endl;
	}
	ofs.close();
}

//write vector to file
template <typename U>
void write2file(vector<U> v, string filename)
{
	ofstream ofs(filename, ofstream::out);

	for (size_t i = 0; i < v.size(); i++)
	{
		ofs << v[i] << ";" << endl;
	}
	ofs.close();
}

//write number to file
template <typename V>
void write2file(V a, string filename)
{
	ofstream ofs(filename, ofstream::out);
    ofs << a << endl;
	ofs.close();
}

//fft complex
template <typename M> 
void fft(M& data, int isign)
{
    int n = data.size()/2, nn = n << 1;
    
    int mmax, m, j, istep;
    double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
    
    j = 1;
    for (int i = 1; i < nn; i += 2) 
    {
        if (j > i){
            swap(data[j - 1], data[i - 1]);
            swap(data[j], data[i]);}
        m=n;
        while (m >= 2 && j > m){
            j -= m;
            m >>= 1;}
        j += m;
    }
    
    mmax = 2;
    while (nn > mmax) 
    {
        istep = mmax << 1;
        theta = isign*(pi_2/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        
        for (m = 1; m < mmax; m += 2){
            for (int i = m; i <= nn; i += istep){
                j = i + mmax;
                tempr = wr*data[j - 1] - wi*data[j];
                tempi = wr*data[j] + wi*data[j - 1];
                data[j - 1] = data[i - 1] - tempr;
                data[j] = data[i] - tempi;
                data[i - 1] += tempr;
                data[i] += tempi;
            }
            wr = (wtemp = wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }
        mmax=istep;
    }   
}

//fft real
template <typename M> 
void fftr(M& data, int isign)
{
    int i, i1, i2, i3, i4, n = data.size();
    double c1 = 0.5, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp;
    double theta = M_PI/double(n >> 1);
    
    if (isign == 1){
        c2 = -0.5;
        fft(data, 1);} 
    else{
        c2 = 0.5;
        theta = -theta;}
    
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    for (i = 1; i < (n >> 2); i++)
    {
        i2 = 1 + (i1 = i + i);
        i4 = 1 + (i3 = n - i1);
        h1r = c1*(data[i1] + data[i3]);
        h1i = c1*(data[i2] - data[i4]);
        h2r = -c2*(data[i2] + data[i4]);
        h2i = c2*(data[i1] -data[i3]);
        data[i1] = h1r + wr*h2r - wi*h2i;
        data[i2] = h1i + wr*h2i + wi*h2r;
        data[i3] = h1r - wr*h2r + wi*h2i;
        data[i4] = -h1i + wr*h2i + wi*h2r;
        wr = (wtemp = wr)*wpr -wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
    }
    
    if (isign == 1){
        data[0] = (h1r = data[0]) + data[1];
        data[1] = h1r - data[1];} 
    else{
        data[0] = c1*((h1r = data[0]) + data[1]);
        data[1] = c1*(h1r - data[1]);
        fft(data, -1);}
   
}

//correlation
template <typename M, typename N> 
void correl(M data1, N data2, vector<double>& ans)
{
    int no2, n = data1.size();
    double tmp;
    vector<double> temp(n);
    
    for (int i = 0; i < n; i++) 
    {
        ans[i] = data1[i];
        temp[i] = data2[i];
    }
    
    fftr(ans, 1);
    fftr(temp, 1);
    no2 = n>>1;
    
    for (int i = 2; i < n; i += 2) 
    {
        tmp = ans[i];
        ans[i] = (ans[i]*temp[i] + ans[i + 1]*temp[i + 1])/no2;
        ans[i + 1] = (ans[i + 1]*temp[i] - tmp*temp[i + 1])/no2;
    }
    
    ans[0] = ans[0]*temp[0]/no2;
    ans[1] = ans[1]*temp[1]/no2;
    
    fftr(ans, -1);
}

//window for spectrum 
double window(int i, int n)
{
    return 1.0 - abs(2.0*i/(n - 1.0) - 1.0);
}

//power spectrum
template <typename T>
void powerspec(T data, vector<double>& ans)
{
    int n = ans.size(), m = n >> 1;
    double w, sumw, fac;
    
    for (int i = 0; i < data.size(); ++i)
    { 
        w = window(i, n);
        ans[i] = w*data[i]; 
        sumw += w*w;
    }
    
    fac = 2./sumw;
    fftr(ans, 1);
    
    ans[0] = 0.5*fac*ans[0]*ans[0];
    for(int j = 0; j < m; ++j) { ans[j] = fac*(ans[2*j]*ans[2*j] + ans[2*j + 1]*ans[2*j + 1]); }
    ans[m] = 0.5*fac*ans[m]*ans[m];
    
    ans.erase(ans.begin() + m, ans.end());
}




#endif 
