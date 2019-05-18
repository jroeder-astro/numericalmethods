#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double func1(double &theta, int &m, int &n);
 
double trap_ext_closed(double &a, double &b, 
                       double (*f)(double &, int &, int &), int &m, int &n,
                       int &N); 

double simp_ext_closed(double &a, double &b, 
                       double (*f)(double &, int &, int &), int &m, int &n,
                       int &N);

double func2(double &x, double &p);

double inf_trafo(double &x, double &p);

double inf_int1(double &a, bool &isinf, 
               double (*f)(double &, double &), double &p, int &N); 

double inf_int2(double &a, bool &isinf, 
               double (*f)(double &, double &), double &p, int &N); 

double some_ext_open(double &a, double &b, 
                     double (*f)(double &, double &), double &p, int &N);

double anth_ext_open(double &a, double &b, 
                     double (*f)(double &, double &), double &p, int &N);


int main(){
  double lower = 0; double upper = M_PI/2;
  int N = 100;
  double err = 0;

  // Closed Quadratures

  int n = 1; int m = 1;

  double integral1 = 0;
  for (N = 100; N < 1e+3; N += 100) {
    integral1 = trap_ext_closed(lower, upper, func1, m, n, N);
    printf("I1 for n = %d, m = %d, N = %d trapezoid: %+3.5f\n", 
           n, m, N, integral1); 
  }

  double integral2 = 0;
  for (N = 100; N < 1e+3; N += 100) {
    integral2 = simp_ext_closed(lower, upper, func1, m, n, N);
    printf("I1 for n = %d, m = %d, N = %d simpson's: %+3.5f\n",
           n, m, N, integral2);
  }

  n = 4; m = 2;

  for (N = 100; N < 1e+3; N += 100) {
    integral1 = trap_ext_closed(lower, upper, func1, m, n, N);
    printf("I1 for n = %d, m = %d, N = %d trapezoid: %+3.5f +- %3.5f\n", 
           n, m, N, integral1, err); 
  }

  for (N = 5; N < 40; N *= 2) {
    integral2 = simp_ext_closed(lower, upper, func1, m, n, N);
    err = fabs(integral2-0.025)/integral2;
    printf("I1 for n = %d, m = %d, N = %d simpson's: %+3.5f +- %3.5f\n",
           n, m, N, integral2, err);
  }


  // Open Quadratures

  lower = 0; bool isinf = true; double p = 0.5; N = 1e+6;
  double int_inf1 = inf_int1(lower, isinf, func2, p, N);
  double int_inf2 = inf_int2(lower, isinf, func2, p, N);
  printf("I2 for p = %+3.3f, open ext N^2: %+3.5f\n", p, int_inf1);
  printf("I2 for p = %+3.3f, open ext N^4: %+3.5f\n", p, int_inf2);

  return 0;
}


double func1(double &theta, int &m, int &n){
  return pow(sin(theta), 2*m-1)*pow(cos(theta), 2*n-1);
}

double func2(double &x, double &p){
  return pow(x, p-1.)/(1+x);
}

double inf_trafo(double &x, double &p){
  return  1./(x*x)* pow(1./x, p-1.)/(1+1./x);
}

double inf_int1(double &a, bool &isinf, 
               double (*f)(double &, double &), double &p, int &N){
  double result = 0; double b_temp = 1;
  if (a == 0) {
    result += some_ext_open(a, b_temp, f, p, N);
  }
  result += some_ext_open(a, b_temp, inf_trafo, p, N);
 
  return result;
} 

double inf_int2(double &a, bool &isinf, 
               double (*f)(double &, double &), double &p, int &N){
  double result = 0; double b_temp = 1;
  if (a == 0) {
    result += anth_ext_open(a, b_temp, f, p, N);
  }
  result += anth_ext_open(a, b_temp, inf_trafo, p, N);
 
  return result;
} 


double trap_ext_closed(double &a, double &b, 
                       double (*f)(double &, int &, int &), int &m, int &n, 
                       int &N){
  if (a == b) return 0;
  double step = (b-a)/(double)N; 
  double result = 1/2*(f(a, m, n)+f(b, m, n));

  double lvl = 0;
  for (int i1 = 1; i1 < N; i1++) {
    lvl = a+(double)i1*step;
    result += f(lvl, m, n);
  }

  return result*step;
}

double simp_ext_closed(double &a, double &b, 
                       double (*f)(double &, int &, int &), int &m, int &n, 
                       int &N){
  if (a == b) return 0;
  double step = (b-a)/(double)N;
  double result = 1/3*(f(a, m, n)+f(b, m, n));
  double lvl1 = 0; double lvl2 = 0;

  for (int i1 = 1; i1 < N; i1+=2) {    //    vvvv  half a step?
    lvl1 = a+(double)i1*step; lvl2 = lvl1-step;
    result += 4./3.*f(lvl1, m, n)+2./3.*f(lvl2, m, n);
  }

  return result*step;
}

double some_ext_open(double &a, double &b, 
                     double (*f)(double &, double &), double &p, int &N){
  if (a == b) return 0;
  double step = (b-a)/(double)N; 
  a += step; b -= step;
  double res = 3./2.*(f(a, p)+f(b, p));

  double lvl = 0;
  for (int i1 = 1; i1 < N-2; i1++) {
    lvl = a+(double)i1*step;
    res += f(lvl, p);
  }

  a -= step; b += step;
  return res*step;
}

double anth_ext_open(double &a, double &b, 
                     double (*f)(double &, double &), double &p, int &N){
  if (a == b) return 0;
  double step = (b-a)/(double)N; 
  a += step; b -= step;
  double res = 27./12.*(f(a, p)+f(b, p));
  a += 2*step; b -= 2*step;
  res += 13./12.*(f(a, p)+f(b, p));
 
  double lvl1 = 0; double lvl2 = 0;

  for (int i1 = 1; i1 < N-3; i1+=2) {
    lvl1 = a+(double)i1*step; lvl2 = lvl1-step;
    res += 4./3.*f(lvl1, p)+2./3.*f(lvl2, p);
  }

  a -= 3*step; b += 3*step;
  return res*step;
}


// test stuff for some_ext_open()

//  printf("result1 seo: %+3.5f\n", res);
//  printf("seo1 f(a, p): %+3.5f\n", f(a, p));
//  printf("seo1 f(b, p): %+3.5f\n", f(b, p));
//  printf("seo1 p = %+3.3f\n", p);
//  printf("seo1 a = %+3.5f\n", a);
//  printf("seo1 b = %+3.5f\n", b);
//  printf("seo1 step = %+3.5f\n", step);

// test stuff for inf_int()

//  printf("result1 ii: %+3.5f\n", result);
//  printf("a = %3.5f\n", a);
//  printf("b_temp = %3.5f\n", b_temp);
 
