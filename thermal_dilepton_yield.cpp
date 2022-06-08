#include <iostream>
#include <cmath>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <ctime>

using namespace std;

//**************************************************************************************************************
//Compile using: g++ -Wall -o thermal_dilepton_yield.out thermal_dilepton_yield.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//************************************************************************************************************
//main function
//************************************************************************************************************

double bessi1(double x)
{
    double ax,ans;
    double y;
    if ((ax=fabs(x)) < 3.75) { 
        y=x/3.75;
        y*=y; 
        ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
            +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))); 
    } else {
        y=3.75/ax; 
        ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
            -y*0.420059e-2)); 
        ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
            +y*(0.163801e-2+y*(-0.1031555e-1+y*ans)))); 
        ans *= (exp(ax)/sqrt(ax));
}
return x < 0.0 ? -ans : ans; 
}


double bessk1(double x)
{
    double bessi1(double x);
    double y,ans; 
    if (x <= 2.0) 
    { 
        y=x*x/4.0;
        ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144 
            +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1 
            +y*(-0.110404e-2+y*(-0.4686e-4)))))));
    } else {
        y=2.0/x;
        ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619 
            +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2 
            +y*(0.325614e-2+y*(-0.68245e-3)))))));
    }
    return ans; 
}

#define pi 3.14159265358979323846
#define COLS 12 // Number of columns in data
double M;
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 1.0;

    double f_0 = 2.25;
    // cout << "f_0? " << endl;
    // cin >> f_0;
    // cout << "f_0 = " << f_0 << endl;
    
    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    fstream file;
    vector < vector <double> > energy_density; // 2d array as a vector of vectors
    vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    int row = 0; // Row counter
    file.open("f0="+ std::to_string(f_0) + ".dat");
    if (file.is_open()) { 
        while (file.good()) { 
        energy_density.push_back(rowVector); // add a new row,
        for (int col=0; col<COLS; col++) {
            file >> std::scientific >> energy_density[row][col]; // fill the row with col elements
        }
        row++; // Keep track of actual row 
        }
    }
    else cout << "Unable to open file" << endl;
    file.close();

    cout << energy_density[102][0] << endl;

    // ofstream myfile1;
    // myfile1.open("thermal_dilepton_yield_xi=1.5_f0=" + std::to_string(f_0) + ".dat");

    // for (int i=1; i<=40; i++)
    // {
    //     M = 0.1000*i;

    //     double tau;
    //     double T;
    //     double mu;
    //     double sum=0.0;
    //     double R=0.0;

    //     for (int j=0; j<42; j++)
    //     {
    //         double hbarc = 0.1973269788;
    //         double alpha = 1.0/137.0;
    //         //double A = pi*(1.0*1.0);
    //         double A = 73.90; //fm^2 for 0-20%
    //         double d_tau = 0.1;
    //         //tau = 2.0+(d_tau*j);
    //         tau = 0.0 + energy_density[j*10][0];
    //         T = 0.0 + energy_density[j*10][1];
    //         mu = 0.0 + energy_density[j*10][2];

    //         //cout << energy_density[j][0] << "\t" << energy_density[j][1] << "\t" << energy_density[j][2] << endl;
            
    //         R = (tau/Q_s)*(d_tau/Q_s)*((5.0*alpha*alpha)/(18.0*pi*pi*pi))*M*(T/Q_s)*bessk1(M/(T/Q_s))*exp((2.0*mu)/T)*A*pow(hbarc,-4.0)*(6.0/5.0)*pow(Q_s,4.0);
    //         sum = sum + R;

    //         // double hbarc = 0.1973269788;
    //         // double alpha = 1.0/137.0;
    //         // //double A = pi*(1.0*1.0);
    //         // double A = 73.90; //fm^2 for 0-20%
    //         // double d_tau = 0.1;
    //         // tau = 1.0+(d_tau*j);
    //         // //tau = 0.5 + energy_density[j*10][0]/(Q_s);
    //         // T = 0.5 + energy_density[j*10][1];
    //         // mu = 0.5 + energy_density[j*10][2];

    //         //cout << tau*hbarc << "\t" << T << "\t" << mu << endl;
            
    //         // R = tau*d_tau*((5.0*alpha*alpha)/(18.0*pi*pi*pi))*M*T*bessk1(M/T)*exp((2.0*mu)/T)*A*pow(hbarc,-2.0)*(6.0/5.0);
    //         // sum = sum + R;            
    //     }
    //     myfile1 << std::scientific << M << "\t" << sum << "\n";
    //     cout << std::scientific << M << "\t" << sum << "\n";
    // }
    // myfile1.close();

    // duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    // std::cout<<"Time: "<< duration <<'\n';
}