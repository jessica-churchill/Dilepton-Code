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

#include <string>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <gsl/gsl_integration.h>

using namespace std;

//**************************************************************************************************************
//Compile using: g++ -Wall -o thermal_dilepton_yield_dNdeta.out thermal_dilepton_yield_dNdeta.cpp -lgsl -lgslcblas
//**************************************************************************************************************


//**************************************************************************************************************
#define pi 3.14159265358979323846
const double alpha = 1.0/137.0;

double T;
double mu;

double I(double M, double Qp, double pp, double pz) //integrand
{
    double E, E1;

    E = sqrt((Qp*Qp)+(M*M));
    E1 = sqrt((pp*pp)+(pz*pz));

    if ( fabs(((2.0*E*E1)-(M*M))/(2.0*Qp*pp)) > 1.0 || ((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)) <= 0.0 || ((Qp*Qp)+(pp*pp)+(pz*pz)-(2.0*E*E1)+(M*M)) < 0.0 )
    {
        return 0;
    }

    else 
    {
        double i = ((alpha*alpha)/(6.0*pow(pi,4.0)))*(M*pp*Qp/(E1*sqrt((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0))))*exp(-(sqrt((pp*pp)+(pz*pz))-mu)/T)*exp(-(sqrt((pp*pp)+(pz*pz))-mu)/T);
        return i;
    }

}
//**************************************************************************************************************


//3d integration routine
//************************************************************************************************************
static double Msav, Qpsav, ppsav;
static double (*nrfunc)(double, double, double, double);

double quad4d(double (*func)(double, double, double, double), double M1, double M2)
{
    double f1(double M, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (10000);

    gsl_function integrand;
    integrand.function = &f1;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, M1, M2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    return result;
}

double f1(double M, void *params) 
{
    double f2(double Qp, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (10000);

    gsl_function integrand;
    integrand.function = &f2;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;
    Msav=M;

    gsl_integration_qag (&integrand, 0.0, 3.0, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double f2(double Qp, void *params) 
{
    double f3(double pp, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (10000);

    gsl_function integrand;
    integrand.function = &f3;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;
    Qpsav=Qp;

    gsl_integration_qag (&integrand, 0.01, 3.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double f3(double pp, void *params) 
{
    double f4(double pz, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (10000);

    gsl_function integrand;
    integrand.function = &f4;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;
    ppsav=pp;

    gsl_integration_qag (&integrand, -3.01, 3.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double f4(double pz, void *params)  
{
    return (*nrfunc)(Msav,Qpsav,ppsav,pz); 
}

//************************************************************************************************************
//main function
//************************************************************************************************************

#define COLS 12 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 2.0;

    //double f_0 = 2.66;

    //ofstream myfile1;
    //myfile1.open("thermal_dilepton_yield_dNdeta_corrected_f0=" + std::to_string(f_0) + ".dat", std::ios_base::app);

    double f_0 = 1.00;
        
    //double f0_array [] = {0.15, 0.30, 0.76, 2.66};
    double limit_array [] = {11, 31, 51, 72, 92};

    for (int f=0; f<5; f++)
    {
        double R = 0.0;
        double sum = 0.0;

        double M1 = 0.0;
        double M2 = 3.0;

        double tau;
        double hbarc = 0.1973269788;

        // double f_0 = f0_array[f];
        // cout << "f0 = "<< f_0 << endl;

        double limit = limit_array[f];

        fstream myfile1;
        myfile1.open("thermal_dilepton_yield_dNdeta_2Qs_f0=" + std::to_string(f_0) + ".dat", std::ios_base::app);

        for (int k=0; k<limit; k++)
        {   
            double d_tau = 0.1;
            //tau = 1.0+(d_tau*((double)k));

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

            tau = energy_density[k*10][0];
            T = energy_density[k*10][1];
            mu = energy_density[k*10][2];
            //cout << tau << "\t" << T << "\t" << mu << endl;
        
            R = (tau/Q_s)*(d_tau/Q_s)*quad4d(I, M1, M2)*pow(hbarc,-2.0)*(6.0/5.0)*pow(Q_s,4.0);
            sum = sum + R;

        }
        myfile1 << std::scientific << tau*hbarc << "\t" << sum << "\n";
        cout << std::scientific << tau*hbarc << "\t" << sum << "\n";
        myfile1.close();
    }    


    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
