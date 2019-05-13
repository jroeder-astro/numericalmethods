#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

double function1(double x);

double function2(double x);

double bisection(double (*f)(double), double upper, double lower, double eps);

bool criterium(double one, double two);

double linear(vector<double> *xvals, vector<double> *yvals, 
              int n, double xbar);
 
main(){
  const int N = 99; 
  vector<double> xi(100);
  vector<double> yi(100);
  double upper =  10;
  double lower = -10;
  double root1 = 0.0;
  double root2 = 0.0;
  double eps  = 1e-7;
  int i1 = 0; int i2 = 0;

  // find roots of given function
  // brackets found by plotting the function
  root1 = bisection(function1, 0, -5., eps);   
  root2 = bisection(function1, 3.5, 7., eps); 
  printf("Roots, precision 1e-7:\n  %+3.8f\n  %+3.8f\n\n", root1, root2);

  // fill up point arrays for interpolation
  double step = (upper-lower)/N;
  for (i1 = 0; i1 < N+1; i1++) xi[i1] = lower + (double)i1 * step;   
  for (i1 = 0; i1 < N+1; i1++) yi[i1] = function1(xi[i1]);   

  // linear interpolation
  // find the points xbar lies between
  printf("f(-5)[intp] = %3.6f\n", linear(&xi,&yi,N+1,-5.));
  printf("f(-5)[true] = %3.6f\n", function1(-5.));
  printf("f(+5)[intp] = %3.6f\n", linear(&xi,&yi,N+1,+5.));
  printf("f(+5)[true] = %3.6f\n", function1(+5.));
 







  return 0;
}


double function1(double x){
  return 3 + 200*x - 30*pow(x, 2.) + 4*pow(x, 3.) - pow(x, 4);
}

double function2(double x){
  return cos(M_PI*x*1/4)+1/2*cos(M_PI*x*3/4)-1/2*cos(M_PI*x*1/2);
}


double neville(double xbar, vector<double> *xvals, 
               vector<double> *yvals){
  int i1 = 0; int i2 = 0;

  

}


double linear(vector<double> *xvals, vector<double> *yvals, 
              int n, double xbar){
  int i1 = 0; int i2 = 0;
  double upx = 0.0;  double upy = 0.0;
  double lowx = 0.0; double lowy = 0.0;

  for (i1 = 0; i1 < n; i1++){
    if ((*xvals)[i1] > xbar) {
      lowx = (*xvals)[i1-1]; lowy = (*yvals)[i1-1];
      upx = (*xvals)[i1];    upy = (*yvals)[i1];
      break;
    }
  }
  
  return ((xbar-upx)*lowy + (lowx-xbar)*upy)/(lowx-upx);
}

double bisection(double (*f)(double), double upper, double lower, double eps){
  // result parameters
  double root      = 1.0;
  double root_prev = 1.0;
  double error     = 5.0;  
  int    count     = 0;

  // function evaluations
  double up  = 0.0;
  double low = 0.0;
  double mid = 0.0;

  // bisection loop
  while (error > eps) {
    root_prev = root;
    low = f(lower);
    up  = f(upper);
    mid = f((lower+upper)/2.);

    if (criterium(low, mid)) upper = (lower + upper) / 2.;    
    else if (criterium(mid, up)) lower = (lower + upper)/2.;

    root = (lower + upper) / 2.;
    error = fabs((root/root_prev)-1);
    
    count += 1;
  } 

  return root;
}

bool criterium(double one, double two){
  if (one*two < 0) return true;
  else return false;
}



