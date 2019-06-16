#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace std;

double u0(double x);

double u0_x(double x);

void ftcs_tstep_atj(int N, vector<vector<double>> &y_n, 
                    vector<double> &y_next, int j, double tau, 
                    double alpha, int tmax_i);

void lxfr_tstep_atj(int N, vector<vector<double>> &y_n, 
                    vector<double> &y_next, int j, double tau, 
                    double alpha, int tmax_i);

void frog_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int tmax_i);

void lxwd_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int tmax_i);


int main() {

  int i1 = 0; int i2 = 0; int i3 = 0; int N = 3; 
  int tsteps = 1000; int xsteps = 1000;

  double v = 1.;                     // factor in wave equation

  double xlow = 0.; double xup = 1.; // range for position
  double tlow = 0.; double tup = 1.; // range for time

  vector<vector<vector<double>>> 
        y(tsteps+2, vector<vector<double>>(xsteps+2, vector<double>(N)));
        // y[i][j][k] <=> index n, index j, and r/s/u
  vector<double> t(tsteps+2); // time array
  vector<double> x(xsteps+2); // position array

  // fill up time array
  for (i1 = 0; i1 <= tsteps+1; i1++) 
    t[i1] = (double)i1 * (tup-tlow)/(double)tsteps;

  // fill up position array and calculate initial values (iv's @ n=0)
  for (i1 = 0; i1 <= xsteps+1; i1++) {
    x[i1] = (double)i1 * (xup-xlow)/(double)xsteps;
    y[0][i1][0] = v * u0_x(x[i1]);
    y[0][i1][1] = 0.;
    y[0][i1][2] = u0(x[i1]);
    // printf("%f, %f, %f, %d\n", v*u0_x(x[i1]), y[0][i1][1], u0(x[i1]), i1);
    // printf("%f, %f, %f, %d\n", y[0][i1][0], y[0][i1][1], y[0][i1][2], i1);
  }

  double alpha = v * (tup-tlow)/tsteps * xsteps/(xup-xlow);
  vector<double> y_t(N);
  double t_unit = (tup-tlow)/tsteps;

  // calculate one time step for every then move to the next time step,
  // again for all x and so on
  for (i1 = 0; i1 <= tsteps; i1++) {
    for (i2 = 0; i2 <= xsteps+1; i2++) {

      // Comment in the one that you want to use 
      // FTCS, Lax-Friedrichs, Lax-Wendroff

      // ftcs_tstep_atj(N, y[i1], y_t, i2, t_unit, alpha, tsteps);
      // lxfr_tstep_atj(N, y[i1], y_t, i2, t_unit, alpha, tsteps);
      lxwd_tstep_atj(N, y[i1], y_t, i2, t_unit, alpha, tsteps);

      for (i3 = 0; i3 < N; i3++) {
        y[i1+1][i2][i3] = y_t[i3];
        // test print
        // printf("y[%d][%d][%d] = %f\n", i1+1, i2, i3, y_t[i3]);
      }
    }
  }

  // output
  for (i3 = 1; i3 <= xsteps; i3++) {
    printf("%f,%f,%f,%d\n", x[i3], y[0][i3][2], y[900][i3][2], i3);
  }

  return 0;
}


double u0(double x) {
  // initial gaussian
  return exp(-pow((x-0.5), 2.) /0.01);
}

double u0_x(double x) {
  // spatial derivative of gaussian, needed for initial conditions
  return (-2./0.01 * (x-0.5))*exp(-pow((x-0.5), 2.) /0.01);
}

void ftcs_tstep_atj(int N, vector<vector<double>> &y_n, 
                    vector<double> &y_next, int j, double tau, 
                    double alpha, int tmax_i) {
  // calculates one time step with FTCS at position j
  // periodic boundary conditions (hopefully)
  int a = j-1; int b = j+1;
  if (j-1 < 1) a = tmax_i;
  if (j+1 > tmax_i+1) b = 1; 

  y_next[0] = y_n[j][0] + alpha/2. * (y_n[b][1] - y_n[a][1]);
  y_next[1] = y_n[j][1] + alpha/2. * (y_n[b][0] - y_n[a][0]);
  y_next[2] = y_n[j][2] + alpha/2. * tau * y_n[j][1];
}

void lxfr_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int tmax_i) {
  // calculates one time step with Lax-Friedrichs at position j
  // periodic boundary conditions (hopefully)
  int a = j-1; int b = j+1;
  if (j-1 < 1) a = tmax_i;
  if (j+1 > tmax_i+1) b = 1; 

  y_next[0] = 1./2.*(y_n[b][0] - y_n[a][0]) + 
              alpha/2. * (y_n[b][1] - y_n[a][1]);
  y_next[1] = 1./2.*(y_n[b][1] - y_n[a][1]) +
              alpha/2. * (y_n[b][0] - y_n[a][0]);
  y_next[2] = y_n[j][2] + alpha/2. * tau * y_n[j][1];
}

void frog_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int tmax_i) {
  // calculates one time step with Leapfrog at position j
  // periodic boundary conditions (hopefully)
  int a = j-1; int b = j+1;
  if (j-1 < 1) a = tmax_i;
  if (j+1 > tmax_i+1) b = 1; 

  // What is the deal with the n-1 point needed?

  // y_next[2] = alpha*alpha*(y_n[b][2] + y_n[a][2]) - y_n[j]
  //            alpha/2. * tau * y_n[j][1];

}

void lxwd_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int tmax_i) {
  // calculates one time step with Lax-Wendroff at position j
  // periodic boundary conditions (hopefully)
  int a = j-1; int b = j+1;
  if (j-1 < 1) a = tmax_i;
  if (j+1 > tmax_i+1) b = 1; 

  y_next[0] = y_n[j][0] + alpha* (1./2.*(y_n[b][1] - y_n[a][1]) + 
              alpha/2. * (y_n[b][0] + y_n[a][0] - 2.* y_n[j][0]));
  y_next[1] = y_n[j][1] + alpha* (1./2.*(y_n[b][0] - y_n[a][0]) + 
              alpha/2. * (y_n[b][1] + y_n[a][1] - 2.* y_n[j][1]));
  y_next[2] = y_n[j][2] + 1./2. * tau * (y_n[j][1] + y_next[1]);
             // Equation (5.15), hopefully s[j][n+1] ^^^^^^^^^
}


