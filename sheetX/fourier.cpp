#include<math.h>
#include<stdlib.h>
#include<stdio.h>
using namespace std;

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

double func(double t);
void four1(float data[], unsigned long nn, int isign);
void print_vec(float vec[], int length);
 
int main(void) {
    int i1 = 0; int i2 = 0; double a = 0.5;
    const int N = 128; const double Del = 1/(double)N;
    float h_n[N]; double sigma = 0.0; int sign = 1;

    for (i1 = 0; i1 < N; i1++) {
        sigma = (double)rand()/(double)RAND_MAX * 2 - 1;
        h_n[i1] = func((double)i1*Del) + a*sigma;
        h_n[i1+1] = (double)i1*Del;
        // printf("%d,%f,%f\n", i1, (double)i1*Del, h_n[i1]);
    }
    
    printf("\n** input vector **\n\n");
    print_vec(h_n, N);

    four1(h_n-1, N, sign);

    printf("\n** output vector **\n\n");
    print_vec(h_n, N);

    printf("\n** exit program **\n\n");

 //   return EXIT_SUCCESS;
}

void print_vec(float vec[], int length) {
    for (int i = 0; i < length; i++) {
        printf("%f\n", vec[i]);
    }
}

double func(double t){
    return -0.5*cos(4*M_PI*t)+0.5*cos(6*M_PI*t)+0.75*cos(18*M_PI*t);
}

void four1(float data[], unsigned long nn, int isign) {
   unsigned long n, mmax, m, j, istep, i;
   double wtemp, wr, wpr, wpi, wi, theta;
   float tempr, tempi;

   n = nn << 1;
   j = 1;

   for (i = 1; i < n; i += 2) {
       if (j > i) {
           SWAP(data[j],data[i]);
           SWAP(data[j+1],data[i+1]);
       }
       m = n >> 1;
       while (m >= 2 && j > m) {
           j -= m;
           m >>= 1;
       }
       j += m;
   }

   mmax = 2;

   while (n > mmax) {
       istep = mmax << 1;
       theta = isign*(6.28318530717959/mmax);
       wtemp = sin(0.5*theta);
       wpr = -2.0*wtemp*wtemp;
       wpi = sin(theta);
       wr = 1.0;
       wi = 0.0;

       for (m = 1; m < mmax; m += 2) {
           for (i = m; i <= n; i += istep) {
               j = i+mmax;
               tempr = wr*data[j] - wi*data[j+1];
               tempi = wr*data[j+1] + wi*data[j];
               data[j]   = data[i]-tempr;
               data[j+1] = data[i+1]-tempi;
               data[i]   += tempr;
               data[i+1] += tempi;
           }
           wr = (wtemp=wr)*wpr - wi*wpi + wr;
           wi = wi*wpr + wtemp*wpi + wr;
       } 
       mmax = istep;
   }
}

