#include<stdio.h> 
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<vector>
using namespace std;

void rk_step(int N, vector<double> *y_chi, double chi, 
             double *y_chi_next, double tau,
             double e, double r1);
void f_times_tau(vector<double> *y_chi, double chi, double *ftau, 
                 double tau, double e, double r1);
double lat(double e, double r1);

main(){

  int N = 2; int i1 = 0; int i2 = 0;
  int max_steps = 2000;
  double tau = 0.01;  

  // remember to put in e and r1 

//  double e = 0.3; double r1 = 7.0;
//  double e = 0.2; double r1 = 6.1;
  double e = 0.22722; double r1 = 5.26;
 

  double y_in[N] = {0.0, 0.0};
  double y_tau[N];
  vector<double> chi(max_steps+1);
  vector<vector<double>> y(max_steps+1, vector<double>(N));
  vector<double> rad(max_steps+1);

  for (i1 = 0; i1 < N; i1++) y[0][i1] = y_in[i1];

  for (i1 = 1; i1 < max_steps; i1++) {
    rk_step(N, &y[i1], chi[i1], y_tau, tau, e, r1);
    for (i2 = 0; i2 < N; i2++) y[i1+1][i2] = y_tau[i2];
    // cout << "b4 chi incr\n";
    chi[i1+1] = chi[i1] + tau;
  } 

  for (i1 = 0; i1 <= max_steps; i1++)
    rad[i1] = lat(e, r1)*5/(1+e*cos(chi[i1]));
  
  // printf("t    phi    r\n");
  for (i1 = 0; i1 <= max_steps; i1++) {
    printf("%5.8lf,%5.8lf,%5.8lf,%5.8lf\n", 
    y[i1][0], y[i1][1], rad[i1], chi[i1]);  
  }

  return 0;
}

double lat(double e, double r1){
  return (1+e)*r1;
}

void rk_step(int N, vector<double> *y_chi, double chi, 
             double *y_chi_next, double tau,
             double e, double r1){
  int i1 = 0;

  double k1[N]; f_times_tau(y_chi, chi, k1, tau, e, r1);

  vector<double> y_tmp(N); 
  for (i1 = 0; i1 < N; i1++) y_tmp[i1] = (*y_chi)[i1] + 0.5*k1[i1];
  double k2[N]; f_times_tau(&y_tmp, chi+0.5*tau, k2, tau, e, r1);
  
  for (i1 = 0; i1 < N; i1++) y_chi_next[i1] = (*y_chi)[i1] + k2[i1];
}

void f_times_tau(vector<double> *y_chi, double chi, double *ftau, 
                 double tau, double e, double r1){
                        //   vvv hard coded black hole mass
  ftau[0] = tau * (lat(e, r1)*5/
            ((lat(e, r1)-2-2*e*cos(chi))*pow(1+e*cos(chi), 2.)) * 
            sqrt((lat(e, r1)-2-2*e)*(lat(e,r1)-2+2*e)/
            (lat(e,r1)-6-2*e*cos(chi))));
  // cout << ftau[0] << endl;
  ftau[1] = tau * sqrt(lat(e,r1)/(lat(e,r1)-6-2*e*cos(chi)));
}

