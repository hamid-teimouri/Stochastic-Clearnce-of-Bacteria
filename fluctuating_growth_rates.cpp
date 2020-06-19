//===============================================================================
// Simulation of Stochastic Clearance of Bacteria (with Fluctuating Growth Rates) Using Gillespie Algorithm
// Dr. Hamid Teimouri, Rice University - Department of Chemistry  
//===============================================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <numeric>
#include <cstdlib>
#include <vector>
#include <valarray>
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <ctime>
#include "ran3.h"
#include "cpu_time.h"
#include <cmath>
#pragma hdrstop
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
using std::vector;
using namespace std;
long int dum;
//==================== difinition of parameters
int Nb=2; // population size
int m=1; //inoculum size
int ind=1; // index for growth rate
int MC=1000000; // MC steps 100000
int n=0; 
int i=0;
double lambda1 = 3.0/60; //bacteria growth rate (per min) 3.0/60
double lambda2 = 0.3/60; //0.3/60;
double gama = (lambda1+lambda2)/2 ; //0.1/60 transition rate between lattices
double delta =gama;
double phi = 0.0; //rate of cell death (per min)
double xm=5;
double x=1.0;
int    Nx=100;
double deltax=xm/double(Nx);
double f1=0.0;
double Nl=0.0;
double t=0.0; 
double t_end=100000000; //100000000
double t_next1=0.0; 
double t_next2=0.0;
double fptl=0.0; 
double sumt=0.0;
double r1=0.0; 
double r2=0.0; 
double h1=0.0; 
double h2=0.0; 
double h3=0.0; 
double h0=0.0; 
double g1=0.0; 
double g2=0.0; 
double g3=0.0; 
double g0=0.0;
//=====================================================================================================
//=========================== MAIN CODE ===============================================================
//=====================================================================================================
int main() {
    dum = -time(NULL);
    ran3(&dum);
    double ctime;
    double ctime1;
    double ctime2;

   ofstream q1("f1_MC.txt");  // extinction probability
   ofstream q2("T1_MC.txt");  // extinction time


//=================================MC simulation============================================
    for (int j = 1; j <= Nx; j++) { // loop over number of initial bacteria
        x = j * deltax;
        phi = x * lambda2;
        Nl = 0; 
        sumt = 0.0; 
        fptl = 0.0; 
            for (int MC_steps = 0; MC_steps <= MC; MC_steps++) {    // loop over MC steps
                n = m;  // inoculum size
                ind=1;  // index for growth rate
                fptl = 0.0; 
                h0 = 0; g0 =0;
                h1 = 0; g1 = 0;
                h2 = 0; g2 = 0;
                h3 = 0; g3 = 0;
                t = 0.0; 
                t_next1 = 0; 
                t_next2 = 0; 

                while (t < t_end) {   // loop over time
                    // calculate propensities for each lattice

                     if(ind==2){ // if lattice
                        g1 = lambda2*n;
                        g2 = phi*n;
                        g3 = gama*n;
                        g0 = g1+g2+g3;
                        t_next2 = ((1 / g0) * (log(1 / ran3(&dum))));  // time to next reaction

                     if (n != 0 && n != Nb) {
                          t = t + t_next2;

                         r2 = ran3(&dum);
                         if (r2 < g1/g0) {
                             n = n + 1;
                         }
                         if (r2 >= g1/g0 && r2 < (g1+g2)/g0) {
                             n = n - 1;
                         }
                         if (r2 >= (g1+g2)/g0) {
                             ind = 1;
                         }
                     }
                         if (n == Nb) {
                         break;
                     }
                     if (n == 0) {
                         fptl = t;
                         Nl = Nl + 1;
                         break;
                     }
                
                     } // end of ind 2

                    if(ind==1){    // if lattice 1
                        h1 = lambda1*n;   // growth
                        h2 = phi*n;      // death
                        h3 = delta*n;     // switch to lattice 2
                        h0 = h1+h2+h3;
                        t_next1 = ((1 / h0) * (log(1 / ran3(&dum))));  // time to next reaction
                    if (n != 0 && n != Nb) {
                        t = t + t_next1;
                    // determine next reaction                    
                            //update time
                        r1 = ran3(&dum);
                            // choose reactions
                        if (r1 < h1/h0) {
                                n = n + 1;   // birth
                            }
                            if (r1 >= h1/h0 && r1 < (h1+h2) / h0) {
                                n = n - 1;   //death
                            }
                            if (r1 >= (h1+h2)/h0) {
                                ind=2;
                            }
                        }
                            if (n == Nb) {
                         break;
                     }
                     if (n == 0) {
                         fptl = t;
                         Nl = Nl + 1;
                         break;
                     }
                       } // end of ind 1

                    } // end of time
                    sumt += fptl;
                   // rf=(N1+N2)
                } // end of mc loop
           
                q1<<x<<" "<<Nl/MC<<endl;
                q2<<x<<" "<<sumt/Nl<<endl;

            } // end of loop over x

//===============================================================================================================================
            ctime2 = cpu_time();
            ctime = ctime2 - ctime1;
            cout << "\n";
            cout << "  Elapsed cpu time for main computation:\n";
            cout << "  " << ctime2 << " seconds.\n";

            q1.close();
            q2.close();
            cout << " end\n " << endl;
            cin.get();

            return 0;
        }