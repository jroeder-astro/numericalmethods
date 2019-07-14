#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<iostream>
using namespace std;

double nonlin(double x, double y);

void gauss_seidel(vector<vector<double>> &y_n, double *y_next,
                  vector<vector<double>> &y_npo, int i, int j, double Del);


main() {

    int i1 = 0; int i2 = 0; int i3 = 0;
    int Nt = 10000; int Ns = 32; // timesteps and spacesteps
    double D = M_PI/(double)Ns; 
    vector<vector<vector<double>>> 
        y(4, vector<vector<double>>(Ns+2, vector<double>(Ns+2)));
    double y_tmp = 0.0;

    // initial values on first three time slices
    // the ghost zones are initialized with the boundary condition
    // here, boundary means phi = phi_A = nonlin(x,y)
    for (i1 = 0; i1 < 4; i1++) {
        for (i2 = 1; i2 < Ns+1; i2++) {
            for (i3 = 1; i3 < Ns+1; i3++) {
                y[i1][0][i3] = nonlin(0, (double)i3*D);
                y[i1][Ns+1][i3] = nonlin(((double)Ns+1.)*D, (double)i3*D);
                y[i1][i2][i3] = 1.0;
                // printf("%f\n", y[i1][i2][i3]);
            }
            y[i1][i2][0] = nonlin((double)i2*D, 0);
            y[i1][i2][Ns+1] = nonlin((double)i2*D, ((double)Ns+1.)*D);
        }
    }

    // do the gauss seidel algorithm initially
    for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 1; i2 < Ns+1; i2++) {
            for (i3 = 1; i3 < Ns+1; i3++) {
                gauss_seidel(y[i1], &y_tmp, y[i1+1], i2, i3, D); 
                y[i1+1][i2][i3] = y_tmp;
                // printf("%f\n", y[i1][i2][i3]);
            }
        }
    }

    // do the gauss seidel between slices 2 and 3 to save memory
    int n = -1; int a = 3;
    for (i1 = 0; i1 < Nt; i1++) {
        for (i2 = 1; i2 < Ns+1; i2++) {
            for (i3 = 1; i3 < Ns+1; i3++) {
                gauss_seidel(y[a], &y_tmp, y[a+n], i2, i3, D); 
                y[a+n][i2][i3] = y_tmp;
                //printf("%f, %d, %d\n", y[a+n][i2][i3], a+n, a);
            }
        }
        n *= -1; 
        if (a == 3) a = 2; 
        else if (a == 2) a = 3;
        //printf("n *= -1, %d\n", i1);
    }

    // output
    for (i1 = 1; i1 < Ns+1; i1++) {
        for (i2 = 1; i2 < Ns+1; i2++) {
            printf("%3.8lf,%3.8lf,%3.8lf\n", 
                  (double)i1*D, (double)i2*D, y[3][i1][i2]);
        }
    }

    return 0;
}


double nonlin(double x, double y){
    return -5.*sin(x+2.*y);
}

void gauss_seidel(vector<vector<double>> &y_n, double *y_next,
                  vector<vector<double>> &y_npo, int i, int j, double Del){
    *y_next = 0.25 * (y_n[i+1][j] + y_n[i][j+1] + y_npo[i-1][j] + 
              y_npo[i][j-1]) - pow((Del/2.),2.) * 
              nonlin((double)i*Del, (double)j*Del);
}

