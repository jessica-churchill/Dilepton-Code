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

//******************************************************************************************************************
//Compile using: g++ -Wall -o nonthermal_dilepton_rate_dRdM2.out nonthermal_dilepton_rate_dRdM2.cpp -lgsl -lgslcblas
//******************************************************************************************************************

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

    double hbarc = 0.1973269788;

    if ( fabs(((2.0*E*E1)-(M*M))/(2.0*Qp*pp)) > 1.0 || ((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)) <= 0.0 || ((Qp*Qp)+(pp*pp)+(pz*pz)-(2.0*E*E1)+(M*M)) < 0.0 )
    {
        return 0;
    }

    else 
    {
        double i = 0.5*((5.0*alpha*alpha)/(72.0*pow(pi,5.0)))*((2.0*pp)/E1)*(Qp/sqrt((4.0*Qp*Qp*pp*pp)-pow(((2.0*E*E1)-(M*M)),2.0)))*f(pz,pp)*f(pz_2,pp_2)*pow(hbarc,-4.0);
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

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, Qp1, Qp2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
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

    gsl_integration_qag (&integrand, 0.01, 3.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);

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

    gsl_integration_qag (&integrand, -3.01, 3.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);

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
    std::clock_t start;
    double duration;
    start = std::clock();

    fstream file;
    vector < vector <double> > array; // 2d array as a vector of vectors
    vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    int row = 0; // Row counter

    //file.open("/Users/JessicaChurchill/Desktop/Final_Code/tracking_dat/2.000000_tracking.dat", ios::in);

    double f_0;
    cout << "f_0? " << endl;
    cin >> f_0;
    cout << "f_0 = " << f_0 << endl;

    double time_step;
    cout << "Time step? " << endl;
    cin >> time_step;
    cout << "Time step = " << time_step << endl;

    file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/f0=" + std::to_string(f_0) + "/" + std::to_string(time_step) + "_tracking.dat", ios::in);
    cout << "/Users/JessicaChurchill/Documents/Research/Final_Code/f0=" + std::to_string(f_0) + "/" + std::to_string(time_step) + "_tracking.dat" << endl;
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


    for (int k=0; k<10201; k++)
    {
        fq[k]=array[k][3];
    }

    double Q_s = 1.0;

    ofstream myfile;
    myfile.open("nonthermal_dilepton_rate_dRdM2_1Qs_f0=" + std::to_string(f_0) + "_t=" + std::to_string(time_step) +".dat");

    for (int i=1; i<=40; i++)
    {
        double Qp1 = 0.01;
        double Qp2 = 3.01;

        M = 0.1000*i;

        double R = quad3d(I, Qp1, Qp2);
        myfile << std::scientific << M << "\t" << R*(6.0/5.0)*pow(Q_s,2.0) << "\n";
        cout << std::scientific << M << "\t" << R*(6.0/5.0)*pow(Q_s,2.0) << "\n";

    }
    myfile.close();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
