#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdlib.h>
using namespace std;

double function(double x, double y);

double bilinear(vector<double> *xvals, vector<double> *yvals, 
                int n, double xbar, double ybar);
 
double linear(vector<double> *xvals, vector<double> *yvals, 
              int n, double xbar);
 

int main(){ 

  double upper =  10;
  double lower = -10;
  int i1 = 0; int i2 = 0;
  double xbar = 5; 
  double ybar = 0;
  //  N+1 grid points in a direction
  const double N = 3;  
  vector<double> result(2);

  //  setting up the grid
  double step = (upper-lower)/N; 
  printf("step = %+2.3f\n",step);
  vector<vector<double>> grid;
  vector<double> row(N+2);

  //  i1 is the y index, i2 is the x index
/*
  for (i1 = 0; i1 < N+1; i1++) {
    grid.push_back(row);
    for (i2 = 0; i2 < N+1; i2++)
      grid[i1][i2] = function(lower+(double)(i2)*step, upper-(double)i1*step);
  }
*/
  for (i1 = 0; i1 < N+1; i1++) {
    grid.push_back(row);
    grid[i1][0] = upper-(double)i1*step;
    for (i2 = 1; i2 < N+2; i2++) grid[i1][i2] = lower+(double)i2*step; 
  }

  //  grid test print  
  for (i1 = 0; i1 < N+1; i1++) {
    printf("| ");
    for (i2 = 0; i2 < N+1; i2++) printf("%+3.3f  ", grid[i1][i2]);
    printf("  |\n");
  }

  return 0;
}


double function(double x, double y){
  return pow(x,2.) + pow(y,2.) + 1;
}

double bilinear(vector<vector<double>> **grid, int n, 
                double xbar, double ybar){
  int i1 = 0;        int i2 = 0;
  double t = 0.0;    double u = 0.0;
  double upx = 0.0;  double upy = 0.0;
  double lowx = 0.0; double lowy = 0.0;

  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 1; i2 < n+1; i2++) {
      if (grid[i1][i2] > xbar) break;
    }
    if (grid[i1][0] < ybar) break;
  }
  
  t = (xbar-grid[0][i2])/(grid[0][i2+1]-grid[0][i2]);
  u = (ybar-grid[i1][0])/(grid[i1][i2+1]-grid[i1][0]);

  return 

}
 
double linear(vector<double> *xvals, vector<double> *yvals, 
              int n, double xbar){
  int i1 = 0; int i2 = 0;
  double upx = 0.0;  double upy = 0.0;
  double lowx = 0.0; double lowy = 0.0;

  for (i1 = 0; i1 < n; i1++) {
    if ((*xvals)[i1] > xbar) {
      lowx = (*xvals)[i1-1]; lowy = (*yvals)[i1-1];
      upx = (*xvals)[i1];    upy = (*yvals)[i1];
      break;
    }
  }
  
  return ((xbar-upx)*lowy + (lowx-xbar)*upy)/(lowx-upx);
}

