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
#include <gsl/gsl_integration.h>

using namespace std;

//**************************************************************************************************************
//Compile using: g++ -Wall -o thermal_dilepton_rate_dRdM2.out thermal_dilepton_rate_dRdM2.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//************************************************************************************************************
//main function
//************************************************************************************************************

#define pi 3.14159265358979323846
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

double M;
double T;
double alpha = 1.0/137.035999084;
double hbarc = .19732698041522198110;

double f_ee(double k, void *params)
{
    double w = sqrt((M*M)+(k*k));
    double f = k*((alpha*alpha)/(6.0*pow(pi, 4.0)*(exp(w/T)-1.0)))*(1.0+(2.0*pow(0.0005109895/M,2.0)))*pow((1.0-(4.0*pow(0.0005109895/M,2.0))),0.5)
                *(((2.0*T)/k)*log(((1.0+exp((-0.5*(w+k))/T))/(1.0+exp((-0.5*(w-k))/T))))+1.0);
    if (4.0*pow(0.0005109895/M,2.0) > 1.0)
    {
        return 0.0;
    }
    else
    {
        return f;
    }
}

double f_mumu(double k, void *params)
{
    double w = sqrt((M*M)+(k*k));
    double f = k*((alpha*alpha)/(6.0*pow(pi, 4.0)*(exp(w/T)-1.0)))*(1.0+(2.0*pow(0.105/M,2.0)))*pow((1.0-(4.0*pow(0.105/M,2.0))),0.5)
                *(((2.0*T)/k)*log(abs((1.0+exp((-0.5*(w+k))/T))/(1.0+exp((-0.5*(w-k))/T))))+1.0);
    if (4.0*pow(0.105/M,2.0) > 1.0)
    {
        return 0.0;
    }
    else
    {
        return f;
    }
}

#define COLS 4 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 1.0;
    
    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    // fstream file;
    // vector < vector <double> > energy_density; // 2d array as a vector of vectors
    // vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    // int row = 0; // Row counter
    // file.open("temperature_chempotential.dat");
    // if (file.is_open()) { 
    //     while (file.good()) { 
    //     energy_density.push_back(rowVector); // add a new row,
    //     for (int col=0; col<COLS; col++) {
    //         file >> std::scientific >> energy_density[row][col]; // fill the row with col elements
    //     }
    //     row++; // Keep track of actual row 
    //     }
    // }
    // else cout << "Unable to open file" << endl;
    // file.close();


    for (int j=1; j<8; j++)
    {
        gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);
        
        gsl_function F_ee;
        F_ee.function = &f_ee;
        
        gsl_function F_mumu;
        F_mumu.function = &f_mumu;

        T = 1.0 + 0.5*j;
        cout << "T = " << T << endl;
        F_ee.params = &T;
        F_mumu.params = &T;

        ofstream myfile1;
        ofstream myfile2;
        ofstream myfile3;
        //myfile1.open("thermal_dilepton_rate_dRdM_ee_and_mumu_" + std::to_string(Q_s) + "Qs_T=" + std::to_string(T) + "GeV.dat");
        myfile2.open("thermal_dilepton_rate_dRdM_ee_T=" + std::to_string(T) + "GeV.dat");
        myfile3.open("thermal_dilepton_rate_dRdM_mumu_T=" + std::to_string(T) + "GeV.dat");
        for (int i=1; i<=30; i++)
        {
            double abs_error = 1.0e-10;  
            double rel_error = 1.0e-10;  
            double result_ee;  
            double error_ee;
            double result_mumu;  
            double error_mumu;

            M = 0.1000*i;
            F_ee.params = &M;
            F_mumu.params = &M;

            // double T = (2.179542e-01)*Q_s;
            // double mu = (-2.898354e-01)*Q_s;
            // double T = 0.2031947690782788;
            // double mu = -0.4669923780341542;
            double mu = 0.0;
            

            gsl_integration_qagiu(&F_ee, 0.0, abs_error, rel_error, 1000, work_ptr, &result_ee, &error_ee);
            //gsl_integration_qag(&F_ee, 0.0, 3.0, abs_error, rel_error, 1000, 1, work_ptr, &result_ee, &error_ee);
            gsl_integration_qagiu(&F_mumu, 0.0, abs_error, rel_error, 1000, work_ptr, &result_mumu, &error_mumu);

            double R = ((alpha*alpha)/(24.0*pi*pi*pi))*M*T*bessk1(M/T)*exp((2.0*mu)/T)*pow(hbarc,-4.0);
            double R_ee = ((2.0*alpha*alpha)/(3.0*pi*pi*pi))*M*M*T*bessk1(M/T)*(1.0+(2.0*pow(0.00051/M,2.0)))*pow((1.0-(4.0*pow(0.00051/M,2.0))), 0.5);
            double R_mumu = ((2.0*alpha*alpha)/(3.0*pi*pi*pi))*M*M*T*bessk1(M/T)*(1.0+(2.0*pow(0.105/M,2.0)))*pow((1.0-(4.0*pow(0.105/M,2.0))), 0.5);
            

            //myfile1 << std::scientific << M << "\t" << 2.0*M*R_ee*pow(Q_s,2.0) << "\t" << 2.0*M*R_mumu*pow(Q_s,2.0) << "\t" << 2.0*pi*M*result_ee*pow(hbarc,-4.0) << "\t" << 2.0*pi*M*result_mumu*pow(hbarc,-4.0) << "\n";
            myfile2 << std::scientific << 2.0*pi*M*result_ee*pow(hbarc,-4.0) << "\n";
            myfile3 << std::scientific << 2.0*pi*M*result_mumu*pow(hbarc,-4.0) << "\n";
            //cout << std::scientific << M << "\t" << M*M*R_ee*pow(Q_s,2.0) << "\t" << M*M*R_mumu*pow(Q_s,2.0) << "\t" << M*M*result_ee*pow(hbarc,-4.0) << "\t" << M*M*result_mumu*pow(hbarc,-4.0) << "\n";
            //cout << "ee cross section: " << "\t" << (1.0+(2.0*pow(0.00051/M,2.0)))*pow((1.0-(4.0*pow(0.00051/M,2.0))), 0.5) << endl;
            //cout << "mumu cross section: " << "\t" << (1.0+(2.0*pow(0.105/M,2.0)))*pow((1.0-(4.0*pow(0.105/M,2.0))), 0.5) << endl;

            cout << M << "\t" << (result_ee+result_mumu)*pow(hbarc,-4.0) << "\t" << R << endl;
        }

        gsl_integration_workspace_free(work_ptr);
        myfile1.close();
        myfile2.close();
        myfile3.close();
        
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
