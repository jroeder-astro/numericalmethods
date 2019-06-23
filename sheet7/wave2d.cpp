#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace std;

/*
  --- Questions ---

  - Leapfrog: how to implement efficiently?
  - Boundary: Only for u or also for r and s?

*/


double u0(double x);

double u0_x(double x, double z);

double u0_z(double x, double z);

void lxwd_tstep_atj(int N, vector<vector<double>> &y_n, 
	            vector<double> &y_next, int j, double tau, 
                    double alpha, int tmax_i);


int main() {

  int i1 = 0; int i2 = 0; int i3 = 0; int N = 4; 
  int tsteps = 3000; int xsteps = 1000; int zsteps = 1000;

  double v = 1.;                     // factor in wave equation

  double xlow = 0.; double xup = 1.; // range for x  position
  double zlow = 0.; double zup = 1.; // range for z  position
  double tlow = 0.; double tup = 1.; // range for time

  vector<vector<vector<vector<double>>>>
              y(tsteps+2, vector<vector<vector<double>>>(xsteps+2, 
              vector<vector<double>>(zsteps+2, vector<double>(N))));
              // y[i][j][k][l] <=> i: time, i: x-pos., j: z-pos., k: r/l/s/u
  vector<double> t(tsteps+2); // time array
  vector<double> x(xsteps+2); // x position array
  vector<double> z(zsteps+2); // z position array

  // fill up t, x & z arrays
  for (i1 = 0; i1 <= tsteps+1; i1++) 
    t[i1] = (double)i1 * (tup-tlow)/(double)tsteps;
  for (i1 = 0; i1 <= xsteps+1; i1++) 
    x[i1] = (double)i1 * (xup-xlow)/(double)xsteps;
  for (i1 = 0; i1 <= zsteps+1; i1++) 
    z[i1] = (double)i1 * (zup-zlow)/(double)zsteps;

  // calculate initial values (iv's @ n=0)
  for (i2 = 0; i2 <= zsteps+1; i2++){
    for (i1 = 0; i1 <= xsteps+1; i1++) {
      y[0][i1][i2][0] = v * u0_x(x[i1], z[i2]);  // r(t=0)
      y[0][i1][i2][1] = v * u0_z(x[i1], z[i2]);  // l(t=0)
      y[0][i1][i2][2] = 0;                       // s(t=0)
      y[0][i1][i2][3] = u0(x[i1], z[i2]);        // u(t=0)
    }
  }

  double alpha_x = v * (tup-tlow)/tsteps * xsteps/(xup-xlow);
  double alpha_z = v * (tup-tlow)/tsteps * zsteps/(zup-zlow);
  vector<double> y_t(N);
  double t_unit = (tup-tlow)/tsteps; 
  double Q = (1.-alpha)/(1.+alpha);

  // calculate one time step for every then move to the next time step,
  // again for all x and so on
  for (i1 = 0; i1 <= tsteps; i1++) {
    for (i2 = 1; i2 <= xsteps; i2++) {

      lxwd_tstep_atj(N, y[i1], y_t, i2, t_unit, alpha, xsteps);

      for (i3 = 0; i3 < N; i3++) {
        y[i1+1][i2][i3] = y_t[i3];
        // test print
        // printf("y[%d][%d][%d] = %f\n", i1+1, i2, i3, y_t[i3]);
      }

      // Outgoing wave boundary conditions
      for (i3 = 0; i3 < N; i3++) {
        y[i1+1][0][i3] = y[i1][1][i3] + Q * (y[i1][0][i3]-y[i1+1][1][i3]);
        y[i1+1][xsteps+1][i3] = y[i1][xsteps][i3] + 
                               Q * (y[i1][xsteps+1][i3]-y[i1+1][xsteps][i3]);
      }
    }
  }

  // output
  for (i3 = 1; i3 <= xsteps; i3++) {
    printf("%f,%f,%f,%d\n", x[i3], y[0][i3][2], y[2000][i3][2], i3);
  }

  return 0;
}


double u0(double x, double z) {
  // initial gaussian
  return exp(- (pow(x-0.5, 2.)+pow(z-0.5, 2.))/(0.01+0.01));
}

double u0_x(double x, double z) {
  // spatial x derivative of gaussian
  return (-2./(0.01+0.01)*(x-0.5)) * u0(x, z);
}

double u0_z(double x, double z) {
  // spatial z derivative of gaussian
  return (-2./(0.01+0.01)*(z-0.5)) * u0(x, z);
}

void lxfr_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int tmax_i) {
  // calculates one time step with Lax-Friedrichs at position j
void lxwd_tstep_atj(int N, vector<vector<double>> &y_n, 
	    vector<double> &y_next, int j, double tau, 
	    double alpha, int xmax_i) {
  // calculates one time step with Lax-Wendroff at position j
  int a = j-1; int b = j+1;
/*
  // periodic boundary conditions (hopefully) 
     // vvvvv these work magically
  if (j-1 < 1) a = xmax_i;
  if (j+1 > xmax_i+1) b = 1; 
*/
  // printf("in normal LXWD\n");
  y_next[0] = y_n[j][0] + alpha* (1./2.*(y_n[b][1] - y_n[a][1]) + 
              alpha/2. * (y_n[b][0] + y_n[a][0] - 2.* y_n[j][0]));
  // printf("y_next[0] ok\n");
  y_next[1] = y_n[j][1] + alpha* (1./2.*(y_n[b][0] - y_n[a][0]) + 
              alpha/2. * (y_n[b][1] + y_n[a][1] - 2.* y_n[j][1]));
  y_next[2] = y_n[j][2] + 1./2. * tau * (y_n[j][1] + y_next[1]);
             // Equation (5.15), hopefully s[j][n+1] ^^^^^^^^^
}

