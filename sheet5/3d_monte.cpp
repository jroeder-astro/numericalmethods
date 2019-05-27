#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<vector>
#include<math.h>
using namespace std;

double func1(double x);

void montecarlo(int N, double(*f)(double), vector<double> &out);


main(){

  int N = 0; int i1  = 0;
  int steps[8] = {10, 20, 50, 100, 200, 500, 2000, 5000};
                 //      x        y        z      < min,max >
  vector<double> outer{1., 4., -3., 4., -1., 1.};

  montecarlo(steps[3], func1, outer);  

  return 0;

}

double func1(double x) {
  return x;
}

void montecarlo(int N, double(*f)(double), vector<double> &out){
  int i1 = 0; int N_irr = 0; int N_reg = 0;
//  double result = 0.0; double error = 0.0;
  double ran = 0.0; srand(2);

  vector<double> xyz(3);
  vector<double> xreg(3);
  vector<double> xres(3);

//  double x = 0; double y = 0; double z = 0;
//  double xreg = 0; double yreg = 0; double zreg = 0;
//  double xres = 0; double yres = 0; double zres = 0;
  double div_reg = 0.0;
  double div_res = 0.0;

  // rho has been set to 1.
  
  // very handwaving (laziness)

  div_reg = (out[1]-out[0])*(out[3]-out[2])*(out[5]-out[4]);

  xreg[0] = 1./2.*(pow(out[1], 2.)-pow(out[0], 2.))*
          (out[3]-out[2])*(out[5]-out[4]);
  cout << xreg[0]/div_reg << endl;

//  xreg[0] /= div_reg;

  xreg[1] = 1./2.*(pow(out[3], 2.)-pow(out[2], 2.))*
          (out[1]-out[0])*(out[5]-out[4]);
//  xreg[1] /= div_reg;

  xreg[2] = 1./2.*(pow(out[5], 2.)-pow(out[4], 2.))*
          (out[1]-out[0])*(out[3]-out[2]);
//  xreg[2] /= div_reg;
 
  // Actual shooting

  for (i1 = 0; i1 < N; i1++) {
    xyz[0] = (double)rand()/(double)RAND_MAX * (out[1]-out[0]) + out[0];
    xyz[1] = (double)rand()/(double)RAND_MAX * (out[3]-out[2]) + out[2];
    xyz[2] = (double)rand()/(double)RAND_MAX * (out[5]-out[4]) + out[4];
    cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl; 
    if (xyz[2]*xyz[2]+pow(sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])-3., 2.)<=1.) 
      N_irr++;
    else if (xyz[2]*xyz[2]+pow(sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])-3., 2.)>1.) 
      N_reg++;
  }

  double ratio = (double)N_irr/(double)N_reg; 
  cout << ratio << endl;
  div_res = ratio * div_reg;
  for (i1 = 0; i1 < 3; i1++){
    xres[i1] = (ratio * xreg[i1])/div_res;
    cout<<"xres "<<i1<<" = "<<xres[i1]<<endl;
  }

//  printf("Full result: %3.5f +- %3.5f\n", result, error);
}

