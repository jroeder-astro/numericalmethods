#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*

The residual vector is not the null vector which means that at some
point I got something wrong. The problem are likely the indices.

Probably the issue: L shares the first row with A, U shares the first column
--> maybe no full decomposition

pot. prob. == potentially problematic

*/

// Functions

void arr_alloc(double **arr, int size) {
  *arr = (double *)malloc(size*sizeof(double));
  if (arr == NULL) {
    printf("malloc not successful");
    exit(0);
  }
}

void mat_alloc(double ***mat, int n) {
  // gets us symmetric matrices
  *mat = (double **)malloc(n*sizeof(double *));
  if (mat == NULL) 
    exit(0);

  for (int i1 = 0; i1 < n; i1++) {
    if (((*mat)[i1] = (double *)malloc(n*sizeof(double))) == NULL) 
      exit(0);
  }
}

void lu_decomposition(double **A, double **U, double **L, int n){
  int i1 = 0;
  int i2 = 0;
  int i3 = 0;
  double sum = 0.0;

  L[0][0] = 1;

  // do I need 1's an all the diagonal?
  for (i1 = 0; i1 < n; i1++) {
    L[i1][i1] = 1;
  }

  U[0][0] = A[0][0];

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

void lyb(double **L, double *y, double *b, int n){
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
}

void uxy(double **U, double *x, double *y, int n){
  int i1 = 0;
  int i2 = 0;
  double sum = 0;
  x[n-1] = y[n-1]/U[n-1][n-1];
// ^^^^^   vvvv   vvvvv  pot. prob.
  for (i1 = n-2; i1 = 0; i1--) {
    sum = 0;    // vvvvvvvv pot. prob.
    for (i2 = i1+1; i2 < n; i2++) {
      sum += U[i1][i2]*x[i2];
    }
    x[i1] = (y[i1]-sum)/U[i1][i1];
  }
}

void print_matrix(double **A, int n){
  for (int i1 = 0; i1 < n; i1++) {
printf("| %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f   %+3.3f |\n", 
      A[i1][0], A[i1][1], A[i1][2], A[i1][3], A[i1][4], A[i1][5], A[i1][6], 
      A[i1][7], A[i1][8], A[i1][9]);
  }
}

double residual(double **A, double *x, double *b, double *R, int n) {
  int i1 = 0; int i2 = 0;
  double sum = 0.0;

  for (i1 = 0; i1 < n; i1++) {
    sum = 0;
    for (i2 = 0; i2 < n; i2++) {
      sum += A[i1][i2]*x[i2];
    }
    R[i1] = sum - b[i1];
  }
}

// Main function

int main() {

  double  eps = 1e-5;
  int      i1 = 0;
  int      i2 = 0;
  const int n = 10;

  // *** Memory allocation ***

  double **A = NULL; double **L = NULL; double **U = NULL;
  mat_alloc(&A, n);   mat_alloc(&L, n);  mat_alloc(&U, n);

  double *b = NULL; double *x = NULL; double *y = NULL; double *R = NULL;
  arr_alloc(&b, n); arr_alloc(&x, n); arr_alloc(&y, n); arr_alloc(&R, n);

  // ***  Matrix operation: LU decomposition *** //
 
  // fill up matrix A
  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 0; i2 < n; i2++) {
      A[i1][i2] = 1/(cos((double)i1+eps)+sin((double)i2+eps));
    }
  }

  printf("\n\n*** MATRIX A ***\n\n"); print_matrix(A, n);

  // fill up vector b
  for (i1 = 0; i1 < n; i1++) {
    b[i1] = 1/((double)i1+1) - eps;
  }
  printf("\n\nb vector: b = \n");
  for (i1 = 0; i1 < n; i1++) {
    printf("|  %+.3f  |\n", b[i1]);
  }

  // do decomposition and solve equations
  lu_decomposition(A, U, L, n);
  printf("\n\n*** MATRIX L ***\n\n"); print_matrix(L, n);
  printf("\n\n*** MATRIX U ***\n\n"); print_matrix(U, n);

  lyb(L, y, b, n);
  uxy(U, x, y, n);  

  printf("\n\nResult vector: x = \n");
  for (i1 = 0; i1 < n; i1++) {
    printf("|  %+.3f  |\n", x[i1]);
  }

  // calculate residual vector
  double res = residual(A, x, b, R, n);
  printf("\n\nResidual vector: R = \n");
  for (i1 = 0; i1 < n; i1++) {
    printf("|  %+.3f  |\n", R[i1]);
  }

  // free allocated memory
  for (i1 = 0; i1 < n; i1++) free(A[i1]);
  free(A);  free(b);  free(x);  free(y);  free(R);
  A = NULL; b = NULL; x = NULL; y = NULL; R = NULL;

  return 0;
}

