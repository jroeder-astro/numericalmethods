#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*

I could not get my head around how to properly loop over
Cramer's rule to keep track of operations.

*/

// Functions

void arr_alloc(double **arr, int n) {
  *arr = (double *)malloc(n*sizeof(double));
  if (arr == NULL) {
    printf("malloc not successful");
    exit(0);
  }
}

void mat_alloc(double ***mat, int n) {
  // gets you a symmetric matrix
  *mat = (double **)malloc(n*sizeof(double *));
  if (mat == NULL) {
    printf("malloc not successful");
    exit(0);
  }

  for (int i1 = 0; i1 < n; i1++) {
    if (((*mat)[i1] = (double *)malloc(n*sizeof(double))) == NULL) {
      printf("malloc not successful");
      exit(0);
    }
  }
}

double det3x3(double **mat){
  double det;
  det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2])-
        mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2])+ 
        mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
  // printf("detA = %+f\n", det);
  return det;
}

void cramer(double **mat, double *x, double *b, double det) {

  // HOW DO I LOOP THIS
  x[0] = ((mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2])*b[0]+
         (-mat[0][1]*mat[2][2]+mat[0][2]*mat[2][1])*b[1]+
          (mat[0][1]*mat[1][2]-mat[0][2]*mat[1][1])*b[2])
          /det;
  x[1] = ((-mat[1][0]*mat[2][2]+mat[1][2]*mat[2][0])*b[0]+
           (mat[0][0]*mat[2][2]-mat[0][2]*mat[2][0])*b[1]+
          (-mat[0][0]*mat[1][2]+mat[0][2]*mat[1][0])*b[2])
          /det;
  x[2] = ((mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0])*b[0]+
         (-mat[0][0]*mat[2][1]+mat[0][1]*mat[2][0])*b[1]+
          (mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0])*b[2])
          /det;

}

double residual(double **A, double *x, double *b, double *R, int n) {
  int i1 = 0; int i2 = 0;
  double sum = 0.0;

  for (i1 = 0; i1 < n; i1++) {
    sum = 0;
    for (i2 = 0; i2 < n; i2++) {
      sum += A[i1][i2]*x[i2];
    }
    R[i1] = sum-b[i1];
  }
}

// Main function

int main() {

  double eps = 1e-5;
  int     i1 = 0;
  int     i2 = 0;

  // *** Memory allocation ***

  const int n = 3;
  double  **A = NULL;
  mat_alloc(&A, n);

  double *b = NULL;
  arr_alloc(&b, n);

  double *x = NULL;
  arr_alloc(&x, n);

  double *R = NULL;
  arr_alloc(&R, n);

  // ***  Matrix operation: Cramer's rule *** //
 
  // fill up matrix A
  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 0; i2 < n; i2++) {
      A[i1][i2] = 1/(cos((double)i1+eps)+sin((double)i2+eps));
    }
  }

  // test print matrix A
  printf("\n\n*** MATRIX A ***\n\n");
  for (i1 = 0; i1 < n; i1++) {
    printf("| %+f  %+f  %+f |\n", A[i1][0], A[i1][1], A[i1][2]);
  }

  // fill up vector b
  for (i1 = 0; i1 < n; i1++) {
    b[i1] = 1/((double)i1+1) - eps; 
  }
  printf("b = (%+f, %+f, %+f)\n", b[0], b[1], b[2]);

  // calculate determinant, do Cramer
  double detA = det3x3(A);
  printf("det(A) = %+f\n", detA); 
  cramer(A, x, b, detA);
  printf("x = (%+f, %+f, %+f)\n", x[0], x[1], x[2]);

  // calculate residual vector
  double res = residual(A, x, b, R, n);
  printf("R = (%+f, %+f, %+f)\n", R[0], R[1], R[2]);

  // free allocated memory
  for (i1 = 0; i1 < n; i1++) free(A[i1]);
  free(A); A = NULL;
  free(b); b = NULL;
  free(x); x = NULL;
  free(R); R = NULL;

  return 0;
}

