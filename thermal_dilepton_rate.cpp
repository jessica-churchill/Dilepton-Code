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
//Compile using: g++ -Wall -o thermal_dilepton_rate.out thermal_dilepton_rate.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//**************************************************************************************************************
#define pi 3.14159265358979323846
const double alpha = 1.0/137.0;

double Qp;
const double M = 3.0; //GeV
double T, mu;

double I(double pp, double pz) //integrand
{
    double a1, b1, c1, a2, b2, c2, E, E1;

    E = sqrt((Qp*Qp)+(M*M));
    E1 = sqrt((pp*pp)+(pz*pz));

    a1 = 4.0*E*E;
    a2 = 4.0*E*E;
    b1 = 0.0;
    b2 = 0.0;
    c1 = (4.0*E*E)-pow(((M*M)-(2.0*Qp*pp)),2.0);
    c2 = (4.0*E*E)-pow(((M*M)+(2.0*Qp*pp)),2.0);

    double hbarc = 0.1973269788;

    if ( fabs(((2.0*E*E1)-(M*M))/(2.0*Qp*pp)) > 1.0 || ((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)) <= 0.0 || ((Qp*Qp)+(pp*pp)+(pz*pz)-(2.0*E*E1)+(M*M)) < 0.0 )
    {
        return 0;
    }

    else 
    {
        double i = ((5.0*alpha*alpha)/(72.0*pow(pi,5.0)))*((2.0*pp)/E1)*(1.0/sqrt((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)))*exp(-E/T)*exp((2.0*mu)/T)*pow(hbarc,-4.0);
        return i;
    }

}
//**************************************************************************************************************



//2d integration routine
//************************************************************************************************************
static double ppsav;
static double (*nrfunc)(double, double);

double quad2d(double (*func)(double, double), double pp1, double pp2)
{
    double f1(double pp, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f1;

    double abs_error = 1.0e-10;
    double rel_error = 1.0e-10;
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, pp1, pp2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    return result;
}

double f1(double pp, void *params) 
{
    double f2(double pz, void *params);
    double pz_min(double);
    double pz_max(double);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f2;

    double abs_error = 1.0e-10;
    double rel_error = 1.0e-10;
    double result;  
    double error;

    ppsav=pp;

    gsl_integration_qag (&integrand, -3.0, 3.0, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);

    return result;
}

double f2(double pz, void *params)  
{
    return (*nrfunc)(ppsav,pz); 
}
//************************************************************************************************************

//main function
//************************************************************************************************************

int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 1.0;
    
    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    // T = 2.179542e-01*Q_s; //1Qs
    // mu = -2.898354e-01*Q_s;

    double f0 = 0.15;
    T = 2.034488e-01;
    mu = -4.869558e-01;

    cout << "Q_perp:         " << "log(d^4R/d^4p):     " << endl;
    
    ofstream myfile;
    myfile.open("thermal_dilepton_rate_f0=" + std::to_string(f0) + ".dat");

    for (int i=1; i<=17; i++)
    {
        Qp = 0.1+(0.20*i);

        double pp1 = 0.0130;
        double pp2 = 3.0130;
        double R = quad2d(I, pp1, pp2);

        cout << std::scientific << Qp << "\t" << R*(6.0/5.0)*pow(Q_s,2.0) << endl;
        myfile << std::scientific << Qp << "\t" << R*(6.0/5.0)*pow(Q_s,2.0)  << endl;
    }

    myfile.close();


    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
