#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<iostream>
using namespace std;

double nonlin(double x, double y);

double analytic(double x, double y);

void gauss_seidel(vector<vector<double>> &y_n, double *y_next,
                  vector<vector<double>> &y_npo, int i, int j, double Del);


main() {

    // NOTE: code is very computationally inefficient at this point
    // will rewrite once it works properly (which it doesn't)

    int i1 = 0; int i2 = 0; int i3 = 0;
    int Nt = 10000; int Ns = 32; // timesteps and spacesteps
    double D = M_PI/(double)Ns; 
    vector<vector<vector<double>>> 
        y(Nt, vector<vector<double>>(Ns+2, vector<double>(Ns+2)));
    double y_tmp = 0.0;
    double Lnorm = 0.0; double Ltemp = 0.0; double error = 1e-6;

    // initial for all timesteps
    // the ghost zones are initialized with the boundary condition
    // here, boundary means phi = phi_A = nonlin(x,y)

    for (i1 = 0; i1 < Nt-1; i1++) {
        for (i2 = 1; i2 < Ns+1; i2++) {
            for (i3 = 1; i3 < Ns+1; i3++) {
                y[i1][0][i3] = analytic(0, (double)i3*D);
                y[i1][Ns+1][i3] = analytic(((double)Ns+1.)*D, (double)i3*D);
                y[i1][i2][i3] = 1.0;
                // printf("%f\n", y[i1][i2][i3]);
            }
            y[i1][i2][0] = analytic((double)i2*D, 0);
            // printf("y[%d][%d][%d] = %f\n", i1, i2, 0, y[i1][i2][0]);
            y[i1][i2][Ns+1] = analytic((double)i2*D, ((double)Ns+1.)*D);
        }
    }

    // ********** everything fine up to this point ********** \\
 
    // in the following there are lines with //*** at the ends,
    // those lines belong to a do loop structure instead of the for loop
    // to implement the error check with the L_inf norm

    // do-loop Gauss-Seidel
    for (i1 = 1; i1 < Nt-1; i1++) {
    // i1 = 1;   //***
    // do {      //***

        // updating the new time level with the temporary array
        for (i2 = 1; i2 < Ns+1; i2++) {
            for (i3 = 1; i3 < Ns+1; i3++) {
                gauss_seidel(y[i1], &y_tmp, y[i1+1], i2, i3, D);
                y[i1+1][i2][i3] = y_tmp;

                Ltemp = y_tmp - nonlin((double)i2*D, (double)i3*D);
                if (Ltemp > Lnorm) Lnorm = Ltemp;
                else continue;
            }
        }
        //i1++;                  //***
        //if (i1>Nt-1) break;    //***

    } //while (Lnorm > error);   //***

    // print results

    for (i3 = 0; i3 < Ns+2; i3++) {
        for (i2 = 0; i2 < Ns+2; i2++) {
            printf("%3.8lf,%3.8lf,%3.8lf\n", 
                  (double)i3*D, (double)i2*D, y[1000][i3][i2]);
        }
    }

    return 0;
}


double nonlin(double x, double y){
    return -5.*sin(x+2.*y);
}

double analytic(double x, double y){
    return sin(x+2.*y);
}

void gauss_seidel(vector<vector<double>> &y_n, double *y_next,
                  vector<vector<double>> &y_npo, int i, int j, double Del){
    *y_next = 0.25 * (y_n[i+1][j] + y_n[i][j+1] + y_npo[i-1][j] + 
              y_npo[i][j-1]) - pow((Del/2.),2.) * nonlin(i*Del, j*Del);
}


















