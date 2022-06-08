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
//Compile using: g++ -Wall -o nonthermal_dilepton_yield.out nonthermal_dilepton_yield.cpp -lgsl -lgslcblas
//**************************************************************************************************************

double fq[10201];
int N = 101;
double fd(int i, int j)
{
  return fq[(i*N)+j];
}

double f(double pz, double pp) //distribution function fq
{
    double f_q;

    double pzmin = 1.000000e-02;
    double dpz = 6.000000e-02;
    double ppmin = 1.000000e-02;
    double dpp = 6.000000e-02;

    int i, j;
    i = (int )((fabs(pz)-pzmin)/dpz);
    j = (int )((fabs(pp)-ppmin)/dpp);

    if ( (fabs(pp) > 6.0) || (fabs(pz) > 6.0) )
    {
        f_q = 0;
    }

    else if (i > N || j > N || i+1 > N || j+1 > N) 
    {
        f_q = 0;
    }

    else if (fabs(pp) < ppmin)
    {
        f_q = fd(i,0);
    }

    else
    {
        f_q = (fd(i,j)+fd(i+1,j+1))/2;
    }

    return f_q;
}

//**************************************************************************************************************
#define pi 3.14159265358979323846
const double hbarc = 0.1973269788;
const double alpha = 1.0/137.0;

double M; //GeV

double I(double Qp, double pp, double pz) //integrand
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

    double pp_2 = sqrt((Qp*Qp)+(pp*pp)+(pz*pz)-(2.0*E*E1)+(M*M));
    double pz_2 = -pz;

    if ( fabs(((2.0*E*E1)-(M*M))/(2.0*Qp*pp)) > 1.0 || ((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)) <= 0.0 || ((Qp*Qp)+(pp*pp)+(pz*pz)-(2.0*E*E1)+(M*M)) < 0.0 )
    {
        return 0;
    }

    else 
    {
        double i = ((alpha*alpha)/(6.0*pow(pi,4.0)))*(pp/E1)*(Qp/sqrt((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)))*f(pz,pp)*f(pz_2,pp_2);
        return i;
    }

}
//**************************************************************************************************************



//3d integration routine
//************************************************************************************************************
static double Qpsav, ppsav;
static double (*nrfunc)(double, double, double);

double quad3d(double (*func)(double, double, double), double Qp1, double Qp2)
{
    double f1(double Qp, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f1;

    double abs_error = 1.0e-10;  
    double rel_error = 1.0e-10;  
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, Qp1, Qp2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);
    return result;
}

double f1(double Qp, void *params) 
{
    double f2(double pp, void *params);
    double pp_min(double);
    double pp_max(double);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f2;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;
    Qpsav=Qp;

    gsl_integration_qag (&integrand, 0.01, 6.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double f2(double pp, void *params) 
{
    double f3(double pz, void *params);
    double pz_min(double);
    double pz_max(double);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f3;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;
    ppsav=pp;

    gsl_integration_qag (&integrand, -6.01, 6.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

    return result;
}

double f3(double pz, void *params)  
{
    return (*nrfunc)(Qpsav,ppsav,pz); 
}

//************************************************************************************************************
//main function
//************************************************************************************************************

#define COLS 4 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start_clock;
    double duration;
    start_clock = std::clock();

    double f_0[] = {2.25, 3.81, 5.75, 6.65, 9.5, 11.1, 6.0, 7.0, 10.25, 11.75};

    double xi[] = {1.0, 1.5, 1.0, 1.0, 1.5, 1.5, 1.0, 1.0, 1.5, 1.5};

    double Q_s[] = {1.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0};

    double Qp1 = 0.01;
    double Qp2 = 6.01;

    double tau;

    int start;
    int stop;
        
    for (int f=0; f<10; f++)
    {
        cout << "f_0 = " << f_0[f] << endl; 
        cout << "xi = " << xi[f] << endl; 
        cout << "Q_s = " << Q_s[f] << endl; 

        ofstream file1;
        file1.open("pre-equil_dilepton_yield_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + "_1fmc.dat");

        if (Q_s[f]==1.0)
        {
            start = 0;
            //stop = 11; // 0.4 fm/c
            stop = 41; // 1.0 fm/c
        }
        else 
        {
            start = 0;
            //stop = 31; // 0.4 fm/c
            stop = 91; // 1.0 fm/c
        }

        for (int a=0; a<20; a++)
        {
            //M = 0.1 + (0.1*a);
            M = 0.1 + (0.2*a);
            double sum = 0.0;
            double R = 0.0;

            for (int k=start; k<stop; k++)
            {   
                double d_tau = 0.1;
                tau = 1.0+(d_tau*k);
                //cout << tau << endl;

                fstream file;
                vector < vector <double> > array; // 2d array as a vector of vectors
                vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
                int row = 0; // Row counter

                //file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);
               //file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/xi=1.5_f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);
                file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/f0_" + std::to_string(f_0[f]) + "_xi_" + std::to_string(xi[f]) + "_np_301_nk_64/" + std::to_string(tau) + "_tracking.dat", ios::in);
                if (file.is_open()) { 

                    while (file.good()) { 
                    array.push_back(rowVector); // add a new row,
                    for (int col=0; col<COLS; col++) {
                        file >> std::scientific >> array[row][col]; // fill the row with col elements
                    }
                    row++; // Keep track of actual row 
                    }   
                }
                else cout << "Unable to open file" << endl;
                file.close();

                for (int l=0; l<10201; l++)
                {
                    fq[l]=array[l][3];
                }

                double A;

                if (f<2)
                {
                    A = 100.58; //fm^2 for RHIC 0-20%
                    //A = 114.5; //fm^2 for RHIC 0-10%
                }

                else if (f>=2 && f<6)
                {
                    A = 124.25; //fm^2 for LHC 0-20%
                    //A = 138.5; //fm^2 for LHC 0-10%
                }
                    
                else
                {
                    A = 127.75; //fm^2 for LHC 0-20% @ 5.02 TeV
                    //A = 140.5; //fm^2 for LHC 0-10% @ 5.02 TeV
                }
        
                R = A*(d_tau/Q_s[f])*(tau/Q_s[f])*quad3d(I, Qp1, Qp2)*pow(hbarc,-2.0);
                sum = sum + R;
            }
        file1 << std::scientific << M << "\t" << sum << "\n";
        cout << std::scientific << M << "\t" << sum << "\n";
        }
    file1.close();
    }
    duration = ( std::clock() - start_clock ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
