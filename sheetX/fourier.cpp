#include<math.h>
#include<stdlib.h>
#include<stdio.h>
using namespace std;



/* 
       *** ISSUE ***
    
    h_i should stay the same as it is filled up with values once and from 
    then on only printed as a reference. However, leaving the "float" data
    type for h_i, it changes in a not specified way once h_n is changed by
    the four1 function. Once I change the data type of h_i to "double" which
    should be fine since to my knowledge it is not subject to further use,
    h_i does indeed stay the same for the whole runtime. However, h_n and
    h_c become "nan" after calling the four1 funtion on h_n. I am unable to 
    pin down what exactly is going wrong so for the time being I am happy 
    with a code that does at least anything.

    Compile and run with:
    g++ -o prog fourier.cpp; ./prog

*/



#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

double func(double t);
void four1(float data[], unsigned long nn, int isign);
 
int main() {
    int i1 = 0; int i2 = 0; double a = 0.5;
    const int N = 128; double Del = 1/(double)N;
    float h_n[2*N]; // sampled data
    float h_i[2*N]; // inital data for later comparison
    float h_c[2*N]; // cleaned data
    double sigma = 0.0;

    // filling input array by sampling given function
    for (i1 = 0; i1 < 2*N; i1 += 2) {
        sigma     = (double)rand()/(double)RAND_MAX * 2. - 1.;
        h_n[i1]   = func((float)i1*Del) + (float)a*(float)sigma; // real part
        h_n[i1+1] = (float)i1*Del;                  // imaginary
        h_i[i1]   = func((double)i1*Del) + a*sigma; // real part
        h_i[i1+1] = (double)i1*Del;                 // imaginary
//        h_i[i1]   = h_n[i1];   // real part
//        h_i[i1+1] = h_n[i1+1]; // imaginary
    }


        printf("\nh_i\t\t h_n\t\t h_c\t\t after initialization\n\n");
        for (i1 = 0; i1 < 10; i1++) 
            printf("%f\t %f\t %f\n", h_i[i1], h_n[i1], h_c[i1]);
        printf("...\t\t ...\t\t ...\t\t\n\n");
        printf("adresses:\nh_i\t\t h_n\t\t h_n-1\t\t\n");
        printf("%p\t %p\t %p\n", &h_i, &h_n, &h_n-1);



    // apply fast fourier transformation
    int sign = 1;
    four1(h_n-1, 2*N, sign);


        printf("\nh_i\t\t h_n\t\t h_c\t\t after transform\n\n");
        for (i1 = 0; i1 < 10; i1++) 
            printf("%f\t %f\t %f\n", h_i[i1], h_n[i1], h_c[i1]);
        printf("...\t\t ...\t\t ...\t\t\n\n");
        printf("adresses:\nh_i\t\t h_n\t\t h_c\t\t\n");
        printf("%p\t %p\t %p\n", &h_i, &h_n, &h_c);



    // remove noise and unusual peaks from transform
    double mean = 0.0;
    for (i1 = 0; i1 < 2*N; i1 += 2) mean += h_n[i1];
    mean /= (double)N;
    double nyq = 1/(2*Del);

    for (i1 = 0; i1 < 2*N; i1 += 2) {
        if (fabs(h_n[i1]) < 3*mean)   h_c[i1] = 0; // noise
        else if (fabs(h_n[i1]) > nyq) h_c[i1] = 0; // peaks
        else                          h_c[i1] = h_n[i1];

        h_c[i1+1] = h_n[i1+1];
    }

    
        printf("\nh_i\t\t h_n\t\t h_c\t\t after cleaning\n\n");
        for (i1 = 0; i1 < 10; i1++) 
            printf("%f\t %f\t %f\n", h_i[i1], h_n[i1], h_c[i1]);
        printf("...\t\t ...\t\t ...\t\t\n\n");
        printf("adresses:\nh_i\t\t h_n\t\t h_c\t\t\n");
        printf("%p\t %p\t %p\n", &h_i, &h_n, &h_c);

    
    
    // inverse fourier transform
    sign = -1; // hopefully correct cleaned FT
    four1(h_c-1, 2*N, sign);
    // is this the correct implementation of 
    // Parseval's theorem?
    for (i1 = 0; i1 < 2*N; i1++) h_c[i1]/=(2*N);

    printf("\n** exit program **\n\n");

    return EXIT_SUCCESS;
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

