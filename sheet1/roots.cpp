#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<vector>
using namespace std;

// function headers

double equiltemp(vector<double> *consts, double T);

double derivative(vector<double> *consts, double T);

double bisection(double (*f)(vector<double> *, double), 
                 vector<double> *consts, double upper, 
                 double lower, double eps);

bool criterium(double one, double two);

double newt_raph(double (*f)(vector<double> *, double), 
                 double (*g)(vector<double> *, double), 
                 vector<double> *consts, double upper, 
                 double lower, double eps);

// main function

main(){

  // constants

  double qlc = 1.64030517 * pow(10, 11);   // Qc³/Lambda
  double q   = 2.7 * pow(10., -10.);       // Q
  double lc  = 1.64603517 * pow(10., -1.); // Lambda/c³
  double urc = 1;                          // U/rho*c²
  double T_p = 0.0;  // proton temperature
  double T_g = 0.0;  // photon temperature
  
  vector<double> constants {q, lc, urc, T_p, T_g};  

  // other parameters

  double root_bisection = 0.0;
  double root_newt_raph = 0.0;

  double lower_bracket = 1.6 * pow(10., 7.);
  double upper_bracket = 2.0 * pow(10., 7.);

  double precision = pow(10., -7.);

  // calculation of i)

  constants[3] = pow(10., 9.);
  constants[4] = pow(10., 7.);

  root_bisection = bisection(equiltemp, &constants, upper_bracket, 
                             lower_bracket, precision);

  root_newt_raph = newt_raph(equiltemp, derivative, &constants, upper_bracket,
                             lower_bracket, precision);

  printf("i) The root is at:\n%+10.10e K (bisection)\n", root_bisection);
  printf("%+10.10e K (Newton-Raphson)\n\n", root_newt_raph); 

  // calculation of ii)

  // change of parameters
  constants[3]  = pow(10., 7.);      // new T_p
  constants[4]  = pow(10., 9.);      // new T_g
  constants[2]  = 8.0 * pow(10., -5.); // new urc
  lower_bracket = pow(10., 7.);      // new lower
  upper_bracket = pow(10., 9.);      // new upper

  root_bisection = bisection(equiltemp, &constants, upper_bracket, 
                             lower_bracket, precision);

  root_newt_raph = newt_raph(equiltemp, derivative, &constants, upper_bracket,
                             lower_bracket, precision);

  printf("ii) The root is at:\n%+10.10e K (bisection)\n", root_bisection);
  printf("%+10.10e K (Newton-Raphson)\n", root_newt_raph); 
 
  return 0;
}

// functions

double equiltemp(vector<double> *consts, double T){
//  return (*consts)[0] * (*consts)[1] * pow(T, 3./2.) *
//         (T - (*consts)[3]) / ( (*consts)[2] - T) - 1;
  return (*consts)[0] * (*consts)[2] * sqrt(T) * (T - (*consts)[4]) - 
         (*consts)[1] * (((*consts)[3]/T) - 1);
}

double derivative(vector<double> *consts, double T){
  return (*consts)[0] * (*consts)[2] * (1/(2*sqrt(T)) * (T - (*consts)[4]) 
         + sqrt(T)) + (*consts)[1] * (*consts)[3] / pow(T, 2.);
}

double bisection(double (*f)(vector<double> *, double), vector<double> *consts,
                 double upper, double lower, double eps){
  // result parameters
  double root      = 1.0;
  double root_prev = 1.0;
  double error     = 5.0;  
  int    count     = 0;

  // function evaluations
  double up  = 0.0;
  double low = 0.0;
  double mid = 0.0;

  // control
  // printf("upper: %5.3e\nlower: %5.3e\nprecn: %5.3e\n", upper, lower, eps);

  // bisection loop
  while (error > eps) {
    root_prev = root;
    low = f(consts, lower);
    up  = f(consts, upper);
    mid = f(consts, (lower+upper)/2.);

    if (criterium(low, mid)) {
      upper = (lower + upper) / 2.;
    }

    else if (criterium(mid, up)) {
      lower = (lower + upper)/2.;
    }

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

double newt_raph(double (*f)(vector<double> *, double), 
                 double (*g)(vector<double> *, double), 
                 vector<double> *consts, double upper, 
                 double lower, double eps){
  // result parameters
  double root      = upper;
  double root_prev = 1.0;
  double error     = 5.0;  
  int    count     = 0;

  // Newton-Raphson iteration loop
  while (error > eps) {
    root_prev = root;
    root = root - f(consts, root)/g(consts, root);
    error = fabs((root/root_prev)-1);
    count += 1;
  }

  printf("count, newt_raph: %d\n", count);

  return root;
}

