#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
using namespace std;


double u0(double x);

double u0_x(double x);

void ftcs_tstep_atj(vector<double> &y_n, 
                    vector<double> &y_next, int j, double tau, 
                    double gamma, int tmax_i);

void dufr_tstep_atj(vector<double> &y_n, vector<double> &y_nmo, 
                    vector<double> &y_next, int j, double tau, 
                    double gamma, int tmax_i);

void btcs_tstep(vector<double> &y_n, vector<double> &y_next, 
                vector<vector<double>> &A, int tsteps);

void lu_decomposition(vector<vector<double>> &A, vector<vector<double>> &U, 
                      vector<vector<double>> &L, int n);

void lyb(vector<vector<double>> &L, vector<double> &y, 
         vector<double> &b, int n);

void uxy(vector<vector<double>> &U, vector<double> &x, 
         vector<double> &y, int n);

void print_matrix(vector<vector<double>> &A, int n);

void residual(vector<vector<double>> &A, vector<double> &x, 
              vector<double> &b, vector<double> &R, int n);


int main() {

  int i1 = 0; int i2 = 0; int i3 = 0; int N = 1; 
  int tsteps = 10000; int xsteps = 100;

  double D = 1.;                     // factor in wave equation

  double xlow = 0.; double xup = 1.; // range for position
  double tlow = 0.; double tup = 1.; // range for time

  vector<vector<double>> y(tsteps+2, vector<double>(xsteps+2));
      // y[i][j] <=> index n, index j. value is u
  vector<double> t(tsteps+2); // time array
  vector<double> x(xsteps+2); // position array

  // fill up position array and calculate initial values (iv's @ n=0)
  for (i1 = 0; i1 <= xsteps+1; i1++) {
    x[i1] = (double)i1 * (xup-xlow)/(double)xsteps;
    y[0][i1] = u0(x[i1]);
  }
  
  // double gamma = 2.* D * (tup-tlow)/tsteps * pow((xsteps/(xup-xlow)), 2.);
  vector<double> y_t(N);           // temporary y for FTCS & DUFR
  vector<double> y_temp(xsteps+2); // temporary y for BTCS
  // double t_unit = (tup-tlow)/tsteps; 
  double t_unit = pow((xup-xlow)/xsteps, 2.);
  double gamma = 2.* D * t_unit * pow((xsteps/(xup-xlow)), 2.);
  double Q = (1.-gamma)/(1.+gamma);

  // fill up time array
  for (i1 = 0; i1 <= tsteps+1; i1++) t[i1] = (double)i1 * t_unit;

  // fill up matrix for BTCS
  vector<vector<double>> btcs_mat(xsteps+2, vector<double>(xsteps+2));
  for (i1 = 0; i1 < xsteps+1; i1++) {
    btcs_mat[i1][i1] = 2.* (1.+gamma);
    btcs_mat[i1][i1+1] = -gamma;
    btcs_mat[i1+1][i1] = -gamma;
    btcs_mat[i1+1][i1+1] = 2.* (1.+gamma); 
  }

  // calculate one time step for every then move to the next time step,
  // again for all x and so on
  for (i1 = 0; i1 <= 1000/*tsteps*/; i1++) {

    btcs_tstep(y[i1], y_temp, btcs_mat, xsteps+2);

    // BTCS goes down MEGA FAST (why?!)
    
    for (i2 = 1; i2 <= xsteps; i2++) {
      // Comment in the one that you want to use 
      // FTCS, Du Fort-Frankel, BTCS

      // ftcs_tstep_atj(y[i1], y_t, i2, t_unit, gamma, xsteps);

      // if (i1==0) 
      //    dufr_tstep_atj(y[i1], y[i1], y_t, i2, t_unit, gamma, xsteps);
      // else 
      //    dufr_tstep_atj(y[i1], y[i1-1], y_t, i2, t_unit, gamma, xsteps);

      // y[i1+1][i2] = y_t[0];
      y[i1+1][i2] = y_temp[i2];
    }

    // Periodic boundary conditions
    // y[i1+1][0] = y[i1+1][xsteps];
    // y[i1+1][xsteps+1] = y[i1+1][1];

    // Dirichlet boundary conditions
    // y[i1+1][0] = 0;
    // y[i1+1][xsteps+1] = 0;

    // von Neumann boundary conditions
    y[i1+1][0] = y[i1+1][1];
    y[i1+1][xsteps+1] = y[i1+1][xsteps];
 }

  // output
  for (i3 = 1; i3 <= xsteps; i3++) {
    printf("%f,%f,%f,%d\n", x[i3], y[0][i3], y[5][i3], i3);
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

void ftcs_tstep_atj(vector<double> &y_n, 
                    vector<double> &y_next, int j, double tau, 
                    double gamma, int tmax_i) {
  // calculates one time step with FTCS at position j

  int a = j-1; int b = j+1;
  y_next[0] = y_n[j] + gamma/2. * (y_n[b] + y_n[a] - 2.*y_n[j]);

  // printf("%lf\n", y_next[0]);
}

void dufr_tstep_atj(vector<double> &y_n, vector<double> &y_nmo,
	    vector<double> &y_next, int j, double tau, 
	    double gamma, int tmax_i) {
  // calculates one time step with Du Fort-Frankel at position j

  int a = j-1; int b = j+1;
  double Q1 = (1.-gamma)/(1.+gamma);
  double Q2 = gamma/(1.+gamma);

  y_next[0] = Q1 * y_nmo[j] + Q2 * (y_n[b] + y_n[a]);
}

void btcs_tstep(vector<double> &y_n, vector<double> &y_next, 
        vector<vector<double>> &A, int tsteps) {
  vector<vector<double>> L(tsteps, vector<double>(tsteps));
  vector<vector<double>> U(tsteps, vector<double>(tsteps));
  vector<double> y_int(tsteps); // y_int = U * x in LUx = b

  lu_decomposition(A, U, L, tsteps);
  lyb(L, y_int, y_n, tsteps);
  uxy(U, y_next, y_int, tsteps);
}

void lu_decomposition(vector<vector<double>> &A, vector<vector<double>> &U, 
                      vector<vector<double>> &L, int n) {
  int i1 = 0; int i2 = 0; int i3 = 0; double sum = 0.0;
  // L[0][0] = 1;
  for (i1 = 0; i1 < n; i1++) L[i1][i1] = 1;
  U[0][0] = A[0][0]; // hopefully correct starting point

  for (i2 = 0; i2 < n; i2++) {

    for (i1 = 0; i1 <= i2; i1++) {
      // loop for beta_ij
      sum = 0.0; //vvvvvvv pot. prob.
      for (i3 = 0; i3 < i1; i3++) sum += L[i1][i3]*U[i3][i2];
      U[i1][i2] = A[i1][i2] - sum;
    }
  
    for (i1 = i2+1; i1 < n; i1++) {
      // loop for alpha_ij
      sum = 0.0; //vvvvvvvv pot prob.
      for (i3 = 0; i3 < i2; i3++) sum += L[i1][i3]*U[i3][i2];
      L[i1][i2] = (A[i1][i2] - sum)/U[i2][i2];
    }
  }
}

void lyb(vector<vector<double>> &L, vector<double> &y, vector<double> &b, 
         int n) {
  int i1 = 0; int i2 = 0; double sum = 0;
  y[0] = b[0]/L[0][0];
  
  for (i1 = 1; i1 < n; i1++) {
    sum = 0;      // vvvvvv pot. prob.
    for (i2 = 0; i2 < i1; i2++) sum += L[i1][i2]*y[i2];
    y[i1] = (b[i1]-sum)/L[i1][i1];
  }

  //printf("\n\n Intermediate result: y = \n");
  //for (i1 = 0; i1 < n; i1++) printf("|  %+.10f  |\n", y[i1]);
}

void uxy(vector<vector<double>> &U, vector<double> &x, vector<double> &y, 
         int n) {
  int i1 = 0; int i2 = 0; double sum = 0;
  x[n-1] = y[n-1]/U[n-1][n-1];
// ^^^^^   vvvv   vvvvv  pot. prob.
  i1 = n-2;
       //         vvvvvvvvvv !!!!!!!!!!!!!!!!!!!!
  for (i1 = n-2; i1 >= 0; i1--) {
    sum = 0;    // vvvvvvvv pot. prob.
    for (i2 = i1+1; i2 < n; i2++) sum += U[i1][i2]*x[i2];
    x[i1] = (y[i1]-sum)/U[i1][i1];
  }

  //printf("\n\n Final result: x = \n");
  //for (i1 = 0; i1 < n; i1++) printf("|  %+.10f  |\n", x[i1]);
}

void print_matrix(vector<vector<double>> &A, int n) {
  for (int i1 = 0; i1 < n; i1++) {
  printf("| %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f |\n", 
      A[i1][0], A[i1][1], A[i1][2], A[i1][3], A[i1][4], A[i1][5], A[i1][6], 
      A[i1][7], A[i1][8], A[i1][9]);
  }
}

void residual(vector<vector<double>> &A, vector<double> &x, 
              vector<double> &b, vector<double> &R, int n) {
  int i1 = 0; int i2 = 0;
  double sum = 0.0;

  //printf("\n\nResidual vector: R = \n\n");
  for (i1 = 0; i1 < n; i1++) {
    sum = 0;
    for (i2 = 0; i2 < n; i2++) sum += A[i1][i2]*x[i2];
    R[i1] = sum - b[i1];
    //printf("|  %+.10f  |\n", R[i1]);
  }
}


