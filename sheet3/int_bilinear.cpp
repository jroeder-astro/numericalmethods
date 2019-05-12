#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdlib.h>
using namespace std;

double function(double x);

double linear(vector<double> *xvals, vector<double> *yvals, 
              int n, double xbar);
 

main(){ 

  double upper =  10;
  double lower = -10;
  int i1 = 0; int i2 = 0;
  double xbar = 5; 
  double ybar = 0;
  //  there will be N+1 grid 
  //  points in a direction
  const double N = 3;  

  //  setting up the grid
  double step = (upper-lower)/N; 
  printf("step = %+2.3f\n",step);
  vector<vector<double>> grid;
  vector<double> row(N+2);

  // grid[i1][0] is the y axis, for each
  // grid[i1][i2] for i2 > 0 the x axis is attached

  for (i1 = 0; i1 < N+1; i1++) {
    grid.push_back(row);
    //row[0] = lower+(double)i1*step;
    grid[i1][0] = upper-(double)i1*step;
   
    for (i2 = 1; i2 < N+2; i2++) {
      //if (i2 == 0) grid[i1][i2] = lower+(double)i1*step;
      /*else*/ grid[i1][i2] = lower+(double)(i2-1)*step;
      
      //if (i2 == 0) row[i2] = lower+(double)i1*step;
      /*else*/// row[i2] = lower+(double)i2*step;
      //printf("%3.4f\n", row[i2]);
    } 
   // grid.push_back(row);
   // printf("grid pushed\n");
  }

  // grid test print  
  for (i1 = 0; i1 < N+1; i1++) {
    printf("| ");
    for (i2 = 0; i2 < N+2; i2++) {
      printf("%+3.3f  ", grid[i1][i2]);
    }
    printf("  |\n");
  }

  return 0;
}


double function(double x, double y){
  return pow(x,2.) + pow(y,2.) + 1;
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

