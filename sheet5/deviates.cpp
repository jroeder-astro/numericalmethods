#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<math.h>
#include<vector>
using namespace std;

// function headers

double distribution(double nu, double T);

double derivative(double nu, double T);

void build_dist(double limit, int num, double min, double max,
                vector<vector<double>> *acc, double T);
 
double bisection(double (*f)(double, double), double T, double upper, 
                 double lower, double eps);

bool criterium(double one, double two);


// main function

main(){

  // parameters

  double nu_min = 0;
  double nu_max = 5e+4;
  double alpha  = 0;
  double temperature    = 6e+3;
  double root_bisection = 0.0;
  double lower_bracket  = 9e+3;
  double upper_bracket  = 11e+3;
  double precision      = pow(10., -7.);
//  vector<double> point(2);
  vector<vector<double>> dist_values;
  int N = 10000; int i1 = 0; int i2 = 0;

  root_bisection = bisection(derivative, temperature, upper_bracket,
                             lower_bracket, precision);
//  cout << "root = " << root_bisection << endl;

  alpha = distribution(root_bisection, temperature); 
//  cout << "max(f(nu) = " << alpha << endl;

  build_dist(alpha, N, nu_min, nu_max, &dist_values, temperature);

//  for (i1 = 0; i1 < dist_values.size(); i1++) 
//    printf("%5.5f, %5.5f\n", dist_values[i1][0], dist_values[i1][1]);

  return 0;
}


// functions

double distribution(double nu, double T){
  return 8.*pow(M_PI, 2.)*pow(nu, 2.)*1./(exp(nu/T)-1.);
}

double derivative(double nu, double T){
  return distribution(nu, T) * (2./nu - 1./(T-T*exp(-nu/T)));
}

void build_dist(double limit, int num, double min, double max,
                vector<vector<double>> *acc, double T){
  srand(2);
  int i1 = 0; double albar = 0.0; // alpha_bar
  double nubar = 0.0; double fnubar = 0.0;
  vector<double> point(2);

  for (i1 = 0; i1 < num; i1++) {
    nubar = (double)rand()/(double)RAND_MAX * (max-min) + min;
   // cout << "nubar = " << nubar << endl;
    fnubar = distribution(nubar, T);
   // cout << "fnubar = " << fnubar << endl;
    albar = (double)rand()/(double)RAND_MAX * limit;
   // cout << "albar = " << albar << endl;
  
    if (albar > fnubar) continue;
    else {
      //point[0] = nubar; point[1] = albar;
      //acc->push_back(point);
      printf("%5.5f, %5.5f\n", nubar, albar);
    }
  }
}

/*
double rand_lim(double limit) {
  srand(2);
  double divisor = (double)RAND_MAX/(limit+1);
  double result = 0.0;
  do {result = (double)rand()/divisor} while (result > limit);
  return result;  
}
*/

double bisection(double (*f)(double, double), double T,
                 double upper, double lower, double eps){
  // result parameters
  double root  = 1.0; double root_prev = 1.0;
  double error = 5.0; int    count     = 0;

  // function evaluations
  double up  = 0.0; double low = 0.0; double mid = 0.0;

  // control
  // printf("upper: %5.3e\nlower: %5.3e\nprecn: %5.3e\n", upper, lower, eps);

  // bisection loop
  while (error > eps) {
    root_prev = root;
    low = f(lower, T);
    up  = f(upper, T);
    mid = f((lower+upper)/2., T);

    if (criterium(low, mid)) upper = (lower + upper) / 2.;
    else if (criterium(mid, up)) lower = (lower + upper)/2.;

    root = (lower + upper) / 2.;
    error = fabs((root/root_prev)-1);
    count += 1;
  } 

  printf("count, bisection: %d\n", count);
  return root;
}

bool criterium(double one, double two){
  if (one*two < 0) return true;
  else return false;
}

