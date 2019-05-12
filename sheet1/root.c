#include<stdlib.h>
#include<stdio.h>
#include<math.h>

// functions

double equiltemp(double *constants, double T){
//  return (*consts)[0] * (*consts)[1] * pow(T, 3./2.) *
//         (T - (*consts)[3]) / ( (*consts)[2] - T) - 1;
  double *consts = (double *) constants;
  return (*consts)[0] * (*consts)[2] * sqrt(T) * (T - (*consts)[4]) - 
         (*consts)[1] * (((*consts)[3]/T) - 1);
}

double bisection(double (*f)(double *, double), double *consts,
                 double upper, double lower, double eps){
  // result parameters
  double root = 1.0;
  double root_prev = 1.0;
  double error = 5.0;  

  // function evaluations
  double up  = 0.0;
  double low = 0.0;
  double mid = 0.0;

  // control
  printf("upper: %5.10lf\nlower: %5.10lf\nprec: %5.10lf\n", upper, lower, eps);

  // bisection loop
  while (error > eps) {
    root_prev = root;
    low = equiltemp(consts, lower);
    up  = equiltemp(consts, upper);
    mid = equiltemp(consts, (lower+upper)/2.);

    if (criterium(low, mid)) {
      upper = (lower + upper) / 2.;
    }

    else if (criterium(mid, up)) {
      lower = (lower + upper)/2.;
    }

    root = (lower + upper) / 2.;
    error = fabs((root/root_prev)-1);
  } 

  return root;
}

bool criterium(double one, double two){
  if (one*two < 0) return 1;
  else return 0;
}

/*
double newt_raph(){

}
*/

// main function

int main(){

  // constants

  double qlc = 1.64030517 * pow(10, 11);   // Qc³/Lambda
  double q   = 2.7 * pow(10., -10.);       // Q
  double lc  = 1.64603517 * pow(10., -1.); // Lambda/c³
  double urc = 1;                          // U/rho*c²
  double T_p = pow(10, 9);  // proton temperature
  double T_g = pow(10, 7);  // photon temperature
  
//  vector<double> constants {q, lc, urc, T_p, T_g};  

  double* constants = malloc(5*sizeof(double));

  // other parameters

  double root_bisection = 0.0;
  double root_newt_raph = 0.0;

  double lower_bracket = 1.6 * pow(10., 7.);
  double upper_bracket = 2.0 * pow(10., 7.);

  double precision = pow(10., -7.);

  // calculation

  root_bisection = bisection(equiltemp, &constants, upper_bracket, 
                             lower_bracket, precision);

  printf("The root / equilibrium temperature is at %+10.10lf.\n", 
         root_bisection);  

  free(constants); 
  constants = NULL;

  return 0;
}

