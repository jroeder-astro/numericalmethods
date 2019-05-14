#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

double function2(double x);

 
int main(){
  int N = 99; int i1 = 0; int i2 = 0;
  double upper = 20; double lower = 0;
  step = (upper-lower)/N;
  
  vector<vector<double>> xy;
  vector<double> pt(2);
  for (i1 = 0; i1 < N+1; i1++) {
    pt[0] = lower+(double)i1*step; pt[1] = function2(pt[0]);
    xy.push_back(pt);
  }

  vector<vector<double>> mat;
  vector<double> row(N+1);
  for (i1 = 0; i1 < N+1; i1++) mat.push_back(row);

  vector<double> bvec(N+1);
  vector<double> ddf(N+1);

  mat[0][0] = 1.0; mat[N][N] = 1.0;
  bvec[0] = 0.0; bvec[N] = 0.0;
  
  for (i1 = 1; i1 < N; i1++) {
    mat[i1][i1-1] = (xy[i1][0]   - xy[i1-1][0]) /6.0;
    mat[i1][i1]   = (xy[i1+1][0] - xy[i1-1][0]) /3.0;
    mat[i1][i1+1] = (xy[i1+1][0] - xy[i1][0])   /6.0;

    bvec[i1] = (xy[i1+1][1]-xy[i1][1])/(xy[i1+1][0]-xy[i1][0])-
               (xy[i1][1]-xy[i1-1][1])/(xy[i1][0]-xy[i1-1][0]);
  }

  // solve mat * ddf = bvec with lu_decomposition, lyb and uxb
  // 


  return 0;
}



double function2(double x){
  return cos(M_PI*x*1/4)+1/2*cos(M_PI*x*3/4)-1/2*cos(M_PI*x*1/2);
}

double cA(){
  // either give xy and an index or only the points between which
  // xbar is located
}

void spline(){
  // implement here:
  // double A = cA();
  // return A*y_j + (1_A)*y_j+1 + 1/6*(pow(A, 3)-A)*pow(xj+1-xj,2)+...
}

void lu_decomposition(vector<vector<double>> &A, vector<vector<double>> &U, 
                      vector<vector<double>> &L, int n){
  int i1 = 0;
  int i2 = 0;
  int i3 = 0;
  double sum = 0.0;

  L[0][0] = 1;

  // do I need 1's an all the diagonal?
  for (i1 = 0; i1 < n; i1++) {
    L[i1][i1] = 1;
  }

  U[0][0] = A[0][0]; // hopefully correct starting point

  for (i2 = 0; i2 < n; i2++) {

    for (i1 = 0; i1 <= i2; i1++) {
      // loop for beta_ij
      sum = 0.0; //vvvvvvv pot. prob.
      for (i3 = 0; i3 < i1; i3++) {
        sum += L[i1][i3]*U[i3][i2];
      }

      U[i1][i2] = A[i1][i2] - sum;
    }
  
    for (i1 = i2+1; i1 < n; i1++) {
      // loop for alpha_ij
      sum = 0.0; //vvvvvvvv pot prob.
      for (i3 = 0; i3 < i2; i3++) {
        sum += L[i1][i3]*U[i3][i2];
      }

      L[i1][i2] = (A[i1][i2] - sum)/U[i2][i2];
    }
  }
}

void lyb(vector<vector<double>> *L, vector<double> *y, 
         vector<double> *b, int n){
  int i1 = 0;
  int i2 = 0;
  double sum = 0;
  y[0] = b[0]/L[0][0];
  
  for (i1 = 1; i1 < n; i1++) {
    sum = 0;      // vvvvvv pot. prob.
    for (i2 = 0; i2 < i1; i2++) {
      sum += L[i1][i2]*y[i2];
    }
    y[i1] = (b[i1]-sum)/L[i1][i1];
  }

  printf("\n\n Intermediate result: y = \n");
  for (i1 = 0; i1 < n; i1++) {
    printf("|  %+.10f  |\n", y[i1]);
  }
}

void uxy(vector<vector<double>> *U, vector<double> *x, 
         vector<double> *y, int n){
  int i1 = 0;
  int i2 = 0;
  double sum = 0;
  x[n-1] = y[n-1]/U[n-1][n-1];
// ^^^^^   vvvv   vvvvv  pot. prob.
  i1 = n-2;
  printf("i1 = %d\n", i1);
     //         vvvvvvvvvv !!!!!!!!!!!!!!!!!!!!
  for (i1 = n-2; i1 >= 0; i1--) {
    sum = 0;    // vvvvvvvv pot. prob.
    for (i2 = i1+1; i2 < n; i2++) {
      sum += U[i1][i2]*x[i2];
    }
    printf("do smth");
    x[i1] = (y[i1]-sum)/U[i1][i1];
  }

  printf("\n\n Final result: x = \n");
  for (i1 = 0; i1 < n; i1++) {
    printf("|  %+.10f  |\n", x[i1]);
  }
}

