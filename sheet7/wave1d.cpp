#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace std;

double u0(double x);

double u0_x(double x);

void ftcs_tstep_atj(int N, vector<vector<double>> &y_n, 
                    vector<double> &y_next, double j, 
                    double tau, double alpha);


int main() {

  int i1 = 0; int i2 = 0; int i3 = 0; int N = 3; 
  int tsteps = 1000; int xsteps = 1000;

  double v = 1.; // factor in wave equation

  double xlow = 0.; double xup = 1.; // range for position
  double tlow = 0.; double tup = 3.; // range for time

  vector<vector<vector<double>>> 
        y(tsteps+1, vector<vector<double>>(xsteps+1, vector<double>(N)));
        // y[i][j][k] <=> index n, index j, and r/s/u
  vector<double> t(tsteps+1); // time array
  vector<double> x(xsteps+1); // position array

  // fill up time array
  for (i1 = 0; i1 <= tsteps; i1++) {
    t[i1] = i1 * (tup-tlow)/(double)tsteps;
    // printf("%f, %d\n", t[i1], i1);
  }
  printf("time filled\n");

  // fill up position array and calculate initial values (iv's @ n=0)
  for (i1 = 0; i1 <= xsteps; i1++) {
    x[i1] = i1 * (xup-xlow)/xsteps;
    // printf("position calc\n");
    y[0][i1][0] = v * u0_x(x[i1]);
    y[0][i1][1] = 0.;
    y[0][i1][2] = u0(x[i1]);
  }
  printf("positions and IV's calculated\n");
 
  double alpha = v * (tup-tlow)/tsteps * xsteps/(xup-xlow);
  vector<double> y_t(N);
  double t_unit = (tup-tlow)/tsteps;

  // calculate one time step for every then move to the next time step,
  // again for all x and so on
  for (i1 = 0; i1 <= tsteps; i1++) {
    for (i2 = 0; i2 <= xsteps; i2++) {
      ftcs_tstep_atj(N, y[i1], y_t, x[i2], t_unit, alpha);
      for (i3 = 0; i3 < N; i3++) y[i1+1][i2][i3] = y_t[i3];
    }
  }


  // output? boundaries?


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
                    vector<double> &y_next, double j, 
                    double tau, double alpha) {
  // calculates one time step with FTCS at position j
  y_next[0] = y_n[j][0] + alpha/2. * (y_n[j+1][1] - y_n[j-1][1]);
  y_next[1] = y_n[j][1] + alpha/2. * (y_n[j+1][0] - y_n[j-1][0]);
  y_next[2] = y_n[j][2] + alpha/2. * tau * y_n[j][1];
}


