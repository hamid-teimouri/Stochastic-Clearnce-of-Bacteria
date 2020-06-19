//===============================================================================
// Simulation of Stochastic Clearance of Bacteria Using Gillespie Algorithm
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
int MC=100000; // MC steps
int n=0;
int i=0;
double x=1.0;
double lambda1 = 3.0/60; //bacteria growth rate (per min) 3.0/60
double lambda2 = 0.3/60; //0.3/60;
double lambda = (lambda1+lambda2)/2 ; //0.1/60 transition rate between lattices
double phi = 0 ; //rate of cell death (per min) x*lambda
double xm=5; 
int Nx=100;
double deltax=xm/Nx; 
double f1=0.0;
double Nl=0.0;
double t=0.0;
double t_end=10000000;
double t_next=0.0;
double fptl=0.0;
double sumt=0.0;
double r2=0.0;
double h1=0.0;
double h2=0.0;
double h0=0.0;
double K0=0.0;
double sum0=0.0;
double T1=0.0;
double T2=0.0;
//=====================================================================================================
//=========================== MAIN CODE ===============================================================
//=====================================================================================================
int main()
{
  dum=-time(NULL);
  ran3(&dum);
  double ctime;
  double ctime1;
  double ctime2;

  ofstream q1("f_MC.txt");  // extinction probability
  ofstream q2("T_MC.txt");  // extinction time

//=================================MC simulation============================================

for(int j=1;j<= Nx;j++){
    x=j*deltax;
    phi=x*lambda;
    Nl=0;
    sumt=0.0;
    fptl=0.0;

 for(int MC_steps=0;MC_steps<= MC; MC_steps++){    // loop over MC steps
  n=m; // n fixed
  h0=0;
  h1=0;
  h2=0;
  t=0.0;
  fptl=0.0;
  t_next=0;

  while (t< t_end) {   // loop over time
     // calculate propensities
    h1 = lambda*n;  // cell growth
    h2 = phi*n;    // cell death
    h0 = h1 + h2; // totall propensity
    t_next = ((1/h0)*(log(1/ran3(&dum))));   // time to next reaction/event
    if( n!=0 && n!=Nb){  // in n is no at the boundaries
      //update time
    t = t + t_next;
    r2=ran3(&dum);
    // choose reactions
    if (r2<h1/h0){
      n=n+1;   // cell growth
    }
    if (r2>h1/h0){
      n=n-1;   //  cell death
    }
    }
    if (n==Nb) { // fixation state
      break;
    }
    if (n==0){  // extinction state
      fptl=t;  // first passage time to hit left bounadry (extinction state)
      Nl=Nl+1;  // number of times it hits left boundary (extinction state)
     break; 
    }
  } // end of time
  sumt+=fptl;
 } // end of mc loop
 
 q1<<x<<" "<<Nl/MC<<endl;
 q2<<x<<" "<<sumt/Nl<<endl;

 //cout<<x<<" "<<Nl/MC<<endl;
 //cout<<x<<" "<<sumt/Nl<<endl;

 } // end of loop over x

//===============================================================================================================================
  ctime2 = cpu_time ();
  ctime = ctime2 - ctime1;
  cout << "\n";
  cout << "  Elapsed cpu time for main computation:\n";
  cout << "  " << ctime2 << " seconds.\n";

  q1.close();
  q2.close();

  cout<<" end\n "<<endl;
  cin.get();
  
  return 0;
}
 
