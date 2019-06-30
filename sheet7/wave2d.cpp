#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
using namespace std;


double u0(double x, double z);

double u0_x(double x, double z);

double u0_z(double x, double z);

void lxwd_tstep_atij(int N, vector<vector<vector<double>>> &y_n, 
	    vector<double> &y_next, int i, int j, double tau, 
	    double alx, double aly, int xmax_i, int zmax_j);


int main() {

  int i1 = 0; int i2 = 0; int i3 = 0; int i4 = 0; int N = 4; 
  int tsteps = 200; int xsteps = 100; int zsteps = 100;

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
  double Q_x = (1.-alpha_x)/(1.+alpha_x);
  double Q_z = (1.-alpha_z)/(1.+alpha_z);

  // calculate one time step for every then move to the next time step,
  // again for all x and y and so on
  for (i1 = 0; i1 <= tsteps; i1++) {
    for (i2 = 1; i2 <= xsteps; i2++) {
      for (i3 = 1; i3 <= zsteps; i3++) {

lxwd_tstep_atij(N, y[i1], y_t, i2, i3, t_unit, alpha_x, alpha_z, xsteps, zsteps);

        for (i4 = 0; i4 < N; i4++) {
          y[i1+1][i2][i3][i4] = y_t[i4];
          // test print
          // printf("y[%d][%d][%d] = %f\n", i1+1, i2, i3, y_t[i3]);
        }
      }
    }

    // z ghost zone outgoing wave boundary conditions
    for (i2 = 0; i2 <= xsteps+1; i2++) {
      for (i3 = 0; i3 < N; i3++) {
        y[i1+1][i2][0][i3] = y[i1][i2][1][i3] + Q_z * (y[i1][i2][0][i3] - 
                                                       y[i1+1][i2][1][i3]);
        y[i1+1][i2][zsteps+1][i3] = y[i1][i2][zsteps][i3] + 
                 Q_z * (y[i1][i2][zsteps+1][i3] - y[i1+1][i2][zsteps][i3]);
      }
    }

    // x ghost zone outgoing wave boundary conditions
    for (i2 = 0; i2 <= zsteps+1; i2++) {
      for (i3 = 0; i3 < N; i3++) {
        y[i1+1][0][i2][i3] = y[i1][1][i2][i3] + Q_x * (y[i1][0][i2][i3] - 
                                                       y[i1+1][1][i2][i3]);
        y[i1+1][xsteps+1][i2][i3] = y[i1][xsteps][i2][i3] + 
                 Q_x * (y[i1][xsteps+1][i2][i3] - y[i1+1][zsteps][i2][i3]);
      }
    }
  }

  // output
//  for (i3 = 1; i3 <= xsteps; i3++) {
 //   printf("%f,%f,%f,%d\n", x[i3], z[i3] y[0][i3][2], y[2000][i3][2], i3);
 // }

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

void lxwd_tstep_atij(int N, vector<vector<vector<double>>> &y_n, 
	    vector<double> &y_next, int i, int j, double tau, 
	    double alx, double aly, int xmax_i, int zmax_j) {

  // calculates one time step with Lax-Wendroff at position i,j
  int a = i-1; int b = i+1; int c = j-1; int d = j+1;

  // r l s u <=> 0 1 2 3 in last y_n index

  // eqs. 7.30, 7.31, 7.32
  y_next[0] = y_n[i][j][0] + alx * (1./2.*(y_n[b][j][2] - y_n[a][j][2]) + 
              alx/2. * (y_n[b][j][0] + y_n[a][j][0] - 2.* y_n[i][j][0]));

  y_next[1] = y_n[i][j][1] + aly * (1./2.*(y_n[i][d][2] - y_n[i][c][2]) + 
              aly/2. * (y_n[i][d][1] + y_n[i][c][1] - 2.* y_n[i][j][1]));

  y_next[2] = y_n[i][j][2] + alx * (1./2.*(y_n[b][j][0] - y_n[a][j][0]) + 
              alx/2. * (y_n[b][j][2] + y_n[a][j][2] - 2.* y_n[i][j][2]))+
              aly * (1./2.*(y_n[i][d][1] - y_n[i][c][1]) + aly/2. * 
              (y_n[i][d][2] + y_n[i][c][2] - 2.* y_n[i][j][2]));
 
  // eq. 7.29
  y_next[3] = y_n[i][j][3] + 1./2. * tau * (y_n[i][j][2] + y_next[2]);
}














