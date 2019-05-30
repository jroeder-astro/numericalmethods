#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<vector>
#include<math.h>
using namespace std;

double func1(double x);

double func2(double y);

void montecarlo(int N, double(*f)(double), double min, double max);


main(){

  int N = 0; int i1  = 0;
  double min = 0.; double max = 1.;
  int steps[8] = {10, 20, 50, 100, 200, 500, 2000, 5000};

  printf("\nTrue Result: %3.5f\n\n", M_PI/4.);

  cout << "\n*** FUNCTION 1 ***\n\n";

  for (i1 = 0; i1 < 8; i1++) {
    montecarlo(steps[i1], func1, min, max);
  }

  cout << "\n\n*** FUNCTION 2 ***\n\n";

  for (i1 = 0; i1 < 8; i1++) {
    montecarlo(steps[i1], func2, min, max);
  }

  return 0;
}


double func1(double x){
  return 1./(1.+pow(x, 2.));
}

double func2(double y){
  return func1(y)/((4.-2.*y)/3.);
}

void montecarlo(int N, double(*f)(double), double min, double max){
  int i1 = 0; int i2 = 0;
  double meanfsq  = 0.0; // <f²> 
  double sqmeanf  = 0.0; // <f>²
  double variance = 0.0; 
  double result   = 0.0; double error = 0.0;
  double ran      = 0.0; srand(2);

  for (i1 = 0; i1 < N; i1++) {
    ran = (double)rand()/(double)RAND_MAX * (max-min)+min;
    meanfsq += pow(f(ran), 2.); 
    sqmeanf += f(ran);
  }

  meanfsq /= (double)N; sqmeanf /= (double)N;
  result   = (max-min) * sqmeanf;
  sqmeanf  = sqmeanf*sqmeanf;
  variance = meanfsq - sqmeanf; 
  error = (max-min)*sqrt(variance/(double)N);

  // printf("Number: %d\nResult: %3.5f\nVariance: %3.5f\n", 
  //        N, result, variance);
  printf("Full result: %3.5f +- %3.5f with %d points\n", result, error, N);
}

