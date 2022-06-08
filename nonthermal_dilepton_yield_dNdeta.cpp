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
//Compile using: g++ -Wall -o nonthermal_dilepton_yield_dNdeta.out nonthermal_dilepton_yield_dNdeta.cpp -lgsl -lgslcblas
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
    double dpz = 3.000000e-02;
    double ppmin = 1.000000e-02;
    double dpp = 3.000000e-02;

    int i, j;
    i = (int )((fabs(pz)-pzmin)/dpz);
    j = (int )((fabs(pp)-ppmin)/dpp);

    if ( (fabs(pp) > 3.0) || (fabs(pz) > 3.0) )
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
const double alpha = 1.0/137.0;

double I(double M, double Qp, double pp, double pz) //integrand
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
        double i = ((alpha*alpha)/(6.0*pow(pi,4.0)))*(M*pp*Qp/(E1*sqrt((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0))))*f(pz,pp)*f(pz_2,pp_2);
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

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f1;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, M1, M2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);
    return result;
}

double f1(double M, void *params) 
{
    double f2(double Qp, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

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

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

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

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

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

#define COLS 4 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 2.0;

    //double f_0 = %(f_0)s;

    // cout << "f_0?" << endl;
    // cin >> f_0 >> endl;

    //ofstream myfile1;
    //myfile1.open("final_nonthermal_dilepton_yield_dNdeta.dat", std::ios_base::app);

    //double f0_array [] = {0.15, 0.30, 0.50, 0.70, 0.76, 0.85, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.66, 3.00};
    //double f0_array [] = {0.15, 0.30, 0.76, 2.66};

    double f_0 = 1.00;

    double limit_array [] = {11, 31, 51, 72, 92};

    for (int f=0; f<5; f++)
    {
        double R = 0.0;
        double sum = 0.0;
        double M1 = 0.0;
        double M2 = 3.0;
        double tau, tau1;
        double hbarc = 0.1973269788;
        // double f_0 = f0_array[f];
        // cout << "f0 = "<< f_0 << endl;
        double limit = limit_array[f];
        fstream myfile1;
        myfile1.open("nonthermal_dilepton_yield_dNdeta_2Qs_f0=" + std::to_string(f_0) + ".dat", std::ios_base::app);
        for (int k=0; k<limit; k++)
        {
            double d_tau = 0.1;
            tau = 1.0+(d_tau*k);

            double d_tau1 = 0.1/Q_s;
            tau1 = (1.0+(d_tau*k))/Q_s;

            fstream file;
            vector < vector <double> > array; // 2d array as a vector of vectors
            vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
            int row = 0; // Row counter

            file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);
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
        
            R = tau1*d_tau1*quad4d(I, M1, M2)*pow(hbarc,-2.0)*(6.0/5.0)*pow(Q_s,4.0);
            sum = sum + R;
        }
        myfile1 << std::scientific << tau*hbarc << "\t" << sum << "\n";
        cout << std::scientific << tau*hbarc << "\t" << sum << "\n";
        myfile1.close();
    }


    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
