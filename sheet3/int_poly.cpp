#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

double function1(double x);

double function2(double x);

double neville(double xbar, vector<double> *xvals, vector<double> *yvals, 
               int ordr, int indx);
 
double bisection(double (*f)(double), double upper, double lower, double eps);

bool criterium(double one, double two);

double linear(vector<double> *xvals, vector<double> *yvals, 
              int n, double xbar);
 
int main(){
  int N = 999; 
  vector<double> xi(N+1); vector<double> yi(N+1);
  double upper =  10; double lower = -10;
  double root1 = 0.0; double root2 = 0.0;
  double eps  = 1e-7; int i1 = 0; int i2 = 0;

  // find roots of given function
  // brackets found by plotting the function
  root1 = bisection(function1, 0, -5., eps);   
  root2 = bisection(function1, 3.5, 7., eps); 
  printf("Roots, precision 1e-7:\n  %+3.8f\n  %+3.8f\n\n", root1, root2);

  // fill up point arrays for interpolation
  double step = (upper-lower)/N;
  for (i1 = 0; i1 < N+1; i1++) xi[i1] = lower + (double)i1 * step;   
  for (i1 = 0; i1 < N+1; i1++) yi[i1] = function1(xi[i1]);   

  vector<double> xvals(4); vector<double> yvals(4);
 
  // Neville's algorithm for first function
  upper = 10; lower = 0; double xbar = 5.;
 
  for (i1 = 0; i1 < N+1; i1++) {
    if (xi[i1] > xbar) {
      for (i2 = 0; i2 < 4; i2++) {
        xvals[i2] = xi[i1-2+i2];            
        // printf("xval: %3.4f\n", xvals[i2]);
        yvals[i2] = function1(xi[i1-1+i2]); 
        // printf("yval: %3.4f\n", yvals[i2]);
      }
      break;
    }
  }

  printf("f1(+5)[intp] = %3.8f\n", neville(xbar, &xvals, &yvals, 3, 0));
  printf("f1(+5)[true] = %3.8f\n", function1(xbar));

  upper = 0; lower = -10;xbar = -5.;
 
  for (i1 = 0; i1 < N+1; i1++) {
    if (xi[i1] > xbar) {
      for (i2 = 0; i2 < 4; i2++) {
        xvals[i2] = xi[i1-2+i2];            
        // printf("xval: %3.4f\n", xvals[i2]);
        yvals[i2] = function1(xi[i1-1+i2]); 
        // printf("yval: %3.4f\n", yvals[i2]);
      }
      break;
    }
  }

  printf("f1(-5)[intp] = %3.8f\n", neville(xbar, &xvals, &yvals, 3, 0));
  printf("f1(-5)[true] = %3.8f\n", function1(xbar));


  /*......................*/

  // Neville's algorithm for second function
  // I could have used the same array as above but decided against since
  // I want to try different N

  N = 999; upper = 20; lower = 0;
  step = (upper-lower)/N;
  vector<double> xj(N+1);  vector<double> yj(N+1);
  //vector<double> xvals(4); vector<double> yvals(4);

  for (i1 = 0; i1 < N+1; i1++) xj[i1] = lower + (double)i1*step;   
  for (i1 = 0; i1 < N+1; i1++) yj[i1] = function2(xj[i1]); 
  xbar = 12.;

  for (i1 = 0; i1 < N+1; i1++) {
    if (xj[i1] > xbar) {
      for (i2 = 0; i2 < 4; i2++) {
        xvals[i2] = xj[i1-2+i2];            
        // printf("xval: %3.4f\n", xvals[i2]);
        yvals[i2] = function2(xj[i1-1+i2]); 
        // printf("yval: %3.4f\n", yvals[i2]);
      }
      break;
    }
  }

  printf("f2(+12)[intp] = %3.8f\n", neville(xbar, &xvals, &yvals, 3, 0));
  printf("f2(+12)[true] = %3.8f\n", function2(xbar));

  return 0;
}


double function1(double x){
  return 3 + 200*x - 30*pow(x, 2.) + 4*pow(x, 3.) - pow(x, 4);
}

double function2(double x){
  return cos(M_PI*x*1/4)+1/2*cos(M_PI*x*3/4)-1/2*cos(M_PI*x*1/2);
}

/*
P(0,1,2,3): gesuchtes Polynom
  ^     ^"ordr"
"indx"
*/

double neville(double xbar, vector<double> *xvals, vector<double> *yvals, 
               int ordr, int indx){
  if (ordr-indx == 1) {                        //   vvv
    return ((xbar - (*xvals)[ordr]) * (*yvals)[indx] +
           ((*xvals)[indx] - xbar) * (*yvals)[ordr]) /
           ((*xvals)[indx] - (*xvals)[ordr]);
  } 
                                                        //                  vvv
  return ((xbar - (*xvals)[ordr]) * neville(xbar, xvals, yvals, ordr-1, indx)+
         ((*xvals)[indx] - xbar) * neville(xbar, xvals, yvals, ordr, indx+1))/
         ((*xvals)[indx] - (*xvals)[ordr]);
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



