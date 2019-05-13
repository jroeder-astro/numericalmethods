#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdlib.h>
using namespace std;

double function(double x, double y);

double bilinear(vector<vector<double>> *grid, int n, 
                double xbar, double ybar, vector<double> *coord);


int main(){ 

  double upper =  10;
  double lower = -10;
  int i1 = 0; int i2 = 0;
  double xbar = 5; 
  double ybar = 0;
  //  N+1 grid points in a direction
  const double N = 999;  

  //  In a physical problem the grid will likely be read in as a table,
  //  here we could just call he function four times in the bilinear(...)
  //  piece of code below, but we give it the grid filled with evaluated
  //  points.

  // setting up the grid
  double step = (upper-lower)/N; 
  printf("step = %+2.3f\n",step);
  vector<vector<double>> grid;
  vector<double> row(N+1);
  vector<double> xy(N+1);

  //  i1 is the y index, i2 is the x index
  for (i1 = 0; i1 < N+1; i1++) {
    xy[i1] = lower+(double)(i1)*step;
    grid.push_back(row);
    for (i2 = 0; i2 < N+1; i2++)
      grid[i1][i2] = function(lower+(double)(i2)*step, upper-(double)i1*step);
  }

  //  results
  printf("f(5,0)[intp] = %+3.5f\n", bilinear(&grid, N+1, xbar, ybar, &xy));
  printf("f(5,0)[true] = %+3.5f\n", function(5, 0));

  return 0;
}


double function(double x, double y){
  return pow(x,2.) + pow(y,2.) + 1;
}

double bilinear(vector<vector<double>> *grid, int n, 
                double xbar, double ybar, vector<double> *coord){
  int i1 = 0; int i2 = 0;
  //  find the point to interpolate in the coordinates
  //  this code is only for symmetric, even grids
  for (i1 = 0; i1 < n; i1++) 
    if ((*coord)[i1] > ybar) break;
  for (i2 = 0; i2 < n; i2++) 
    if ((*coord)[i2] > xbar) break;
 
  double t = (xbar-(*coord)[i2-1])/((*coord)[i2]-(*coord)[i2-1]);
  double u = (ybar-(*coord)[i1-1])/((*coord)[i1]-(*coord)[i1-1]);

  return (1-u)*(1-t)*(*grid)[i1-1][i2-1] + t*(1-u)*(*grid)[i1][i2-1] + 
         t*u*(*grid)[i1][i2] + (1-t)*u*(*grid)[i1-1][i2];
}

/*
//  alternative grid filling
  for (i1 = 0; i1 < N+1; i1++) {
    grid.push_back(row);
    grid[i1][0] = upper-(double)i1*step;
    for (i2 = 1; i2 < N+2; i2++) grid[i1][i2] = lower+(double)i2*step; 
  }
*/

/*
//  grid test print  
  for (i1 = 0; i1 < N+1; i1++) {
    printf("| ");
    for (i2 = 0; i2 < N+1; i2++) printf("%+3.3f  ", grid[i1][i2]);
    printf("  |\n");
  }
*/

