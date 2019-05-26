#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<vector>
#include<math.h>

double func1(double x);

void montecarlo(int N, double(*f)(double), double min, double max);


main(){

  int N = 0; int i1  = 0;
  double min = 0.; double max = 1.;
  int steps[8] = {10, 20, 50, 100, 200, 500, 2000, 5000};

  printf("\nTrue Result: %3.5f\n\n", M_PI/4.);

  for (i1 = 0; i1 < 8; i1++) {
    montecarlo(steps[i1], func1, min, max);
  }

  return 0;
}


double func1(double x){
  return 1./(1.+pow(x, 2.));
}

void montecarlo(int N, double(*f)(double), double min, double max){
  int i1 = 0; int i2 = 0;
  double meanfsq  = 0.0; // <f²> 
  double sqmeanf  = 0.0; // <f>²
  double variance = 0.0; double result   = 0.0;
  double ran     = 0.0; srand(2);

  for (i1 = 0; i1 < N; i1++) {
    ran = (double)rand()/(double)RAND_MAX * (max-min)+min;
    meanfsq += pow(f(ran), 2.);
    sqmeanf += f(ran);
  }
  
  meanfsq /= (double)N;
  sqmeanf /= (double)N;
  sqmeanf = sqmeanf*sqmeanf;

  variance = meanfsq - sqmeanf;
  result   = meanfsq + sqmeanf;

  printf("Result: %3.5f\nVariance: %3.5f\nNumber: %d\n", 
         result, variance, N);
}

