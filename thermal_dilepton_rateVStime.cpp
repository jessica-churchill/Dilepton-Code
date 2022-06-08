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

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <gsl/gsl_integration.h>

using namespace std;

//**************************************************************************************************************
//Compile using: g++ -Wall -o thermal_dilepton_rateVStime.out thermal_dilepton_rateVStime.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//**************************************************************************************************************
#define pi 3.14159265358979323846
const double alpha = 1.0/137.0;

const double Qp = 1.0;
const double M = 1.0; //GeV
const double hbarc = 0.1973269788;

double T, mu;

//main function
//************************************************************************************************************
#define COLS 8 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 1.0;
    
    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    double f_0;
    cout << "f_0? " << endl;
    cin >> f_0;
    cout << "f_0 = " << f_0 << endl;

    fstream file;
    vector < vector <double> > energy_density; // 2d array as a vector of vectors
    vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    int row = 0; // Row counter
    file.open("temperature_chempotential_f0="+ std::to_string(f_0) + ".dat");
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
    
    ofstream myfile;
    myfile.open("thermal_dilepton_rateVStime_f0=" + std::to_string(f_0) + ".dat");

    for (int i=0; i<246; i++)
    {
        double tau = energy_density[i*10][0];
        T = energy_density[i*10][1]*Q_s;
        mu = energy_density[i*10][2]*Q_s;

        double E = sqrt((Qp*Qp)+(M*M));
        double R = ((5.0*alpha*alpha)/(72.0*pi*pi*pi*pi))*exp(-E/T)*exp((2.0*mu)/T)*pow(hbarc,-4.0)*(6.0/5.0);

        
        myfile << std::scientific << tau/Q_s << "\t" << R*pow(Q_s,2.0) << "\n";
        cout << std::scientific << tau/Q_s << "\t" << R*pow(Q_s,2.0) << "\n";
    }

    myfile.close();


    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
