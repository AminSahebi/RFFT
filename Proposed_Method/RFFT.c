/************************************************************************************************
 *  Real-Valued FFT										*
 *  Affiliation	: University of Siena & Institute for High Performance Computing and Networking *                                                                          
 *  AUTHOR(s)	: Amin Sahebi, Lorenzo Verdoscia, Roberto Giorgi - UNISI , CNR			*
 *  VERSION	: 0.1.1										*
 *  DATE	: 01/05/2020 (dd/mm/yy)								*                                                 
 *												*
 ************************************************************************************************/
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "time.h"
struct timespec time_start, time_end;

/*******************************************************************************
 * Printing
 ******************************************************************************/
#define logmsg(msg, ...) {\
	    struct timespec _t1;\
	    time_t _secs,_min;\
	    clock_gettime(CLOCK_REALTIME_COARSE,&_t1);\
	    _secs=_t1.tv_sec%60;\
	    _min=(_t1.tv_sec/60)%60;\
	    fprintf(stderr,"[%02d:%02d.%06d]: " msg "\n", (int)_min, (int)_secs, (int)(_t1.tv_nsec/1000), ##__VA_ARGS__);\
}

/*******************************************************************************
 * The Global Variables Definition
 ******************************************************************************/
typedef unsigned char byte;
const float w = 0.707106781187;
//#define DATA double
#define DATA double

DATA t1 = 0;
DATA t2 = 0;

DATA *my, *out;

struct timeval tim;
int Num_Iterations = 100;

int FFT_N; // Number of points in FFT

DATA *X_r, *X_i, *temp_r, *temp_i;

/*******************************************************************************
 * Bit-reversal function to find even and odd numbers
 ******************************************************************************/
unsigned bit_reverse(register int x) {
    x = (x & 0xaaaaaaaa) >> 1 | (x & 0x55555555) << 1;
    x = (x & 0xcccccccc) >> 2 | (x & 0x33333333) << 2;
    x = (x & 0xf0f0f0f0) >> 4 | (x & 0x0f0f0f0f) << 4; 
    x = (x & 0xff00ff00) >> 8 | (x & 0x00ff00ff) << 8;

    return ((x >> 16 ) | (x << 16 ));
}

/*******************************************************************************
 * Real Valued Fast Fourier Transform Body Function
 ******************************************************************************/
void fft(DATA my[], int n) { 
    DATA a_r, a_i;
    unsigned log_n = log2(n);

    for(int i = 0; i < n; i++) {
       unsigned j = bit_reverse(i) >> (sizeof(i) * 8 - log_n); 	  

       	    if (i > j) {
            DATA tmp_r = my[i];
            my[i] = my[j];
            my[j] = tmp_r;
        }
    }

    memset(X_i, 0, FFT_N * sizeof(DATA));
    memset(temp_i, 0, FFT_N * sizeof(DATA));
    memcpy(temp_r, my, FFT_N * sizeof(DATA));

    for (int n = 0; n < FFT_N; n += 2) {
        X_r[n] = temp_r[n] + temp_r[n + 1];
        X_r[n + 1] = temp_r[n] - temp_r[n + 1];
    }

    memcpy(temp_r, X_r, FFT_N * sizeof(DATA));
    memcpy(temp_i, X_i, FFT_N * sizeof(DATA));

    for (int n = 0; n < FFT_N / 2; n++) {
        X_r[n] = temp_r[n] + temp_r[n + FFT_N / 2];
        X_i[n] = temp_i[n] + temp_i[n + FFT_N / 2];
        X_r[n + FFT_N / 2] = temp_r[n] - temp_r[n + FFT_N / 2];
        X_i[n + FFT_N / 2] = temp_i[n] - temp_i[n + FFT_N / 2];
    }
}

/*******************************************************************************
 * Prepare the input and other necessary elements
 ******************************************************************************/
void prepare() {
    my = (DATA*)malloc(FFT_N * sizeof(DATA));
    out = (DATA*)malloc(FFT_N * sizeof(DATA));
    X_r = (DATA*)malloc(FFT_N * sizeof(DATA));
    X_i = (DATA*)malloc(FFT_N * sizeof(DATA));
    temp_r = (DATA*)malloc(FFT_N * sizeof(DATA));
    temp_i = (DATA*)malloc(FFT_N * sizeof(DATA));
    
    for (int i = 0; i < FFT_N; i++) {
        my[i] = (DATA)(1.0 * i) + (DATA)1;	
    }
// srand(time(NULL));

// Loop to generate random numbers in a specified range
//for (int i = 0; i < FFT_N; i++) {
    // Generate a random number between 0 and 1, and scale it to your desired range
  //  double random_value = (double)rand() / RAND_MAX;
   // my[i] = (DATA)(random_value * (100) + 1);
//}
}

/*******************************************************************************
 * Compute the FFT Algorithm, use iterations to normalize the execution time
 ******************************************************************************/
void compute() {
    for (int i = 0; i < Num_Iterations; i++) {
        fft(my, FFT_N);
    }
}


/*******************************************************************************
 * Report the Results
 ******************************************************************************/
void report() {
    DATA sum_input = 0;
    DATA sum_output = 0;

    for (int i = 0; i < FFT_N; i++) {
        sum_input += (my[i] * my[i]);
        sum_output += (X_r[i] * X_r[i] + X_i[i] * X_i[i]);
    }

    // Scale the sum of squares of the output by 1/N
    sum_output /= (DATA)4.0;

    // Calculate the error as the absolute difference
    DATA Err = fabs(sum_input - sum_output);

    logmsg("Using float Variables...");
    logmsg("the abs of inputs are %.9lf", sum_input);
    logmsg("the abs of out puts are %.9lf", sum_output);
    logmsg("The Err for is %.9lf", Err);

    if (Err < 0.0000000000000001) {
        logmsg("****Success****");
    } else {
        logmsg("****Fail****");
    }

    logmsg("Total elapsed time %8ld nsec", time_diff_nsec(&time_end, &time_start) / Num_Iterations);
}

/*****************************************************************************
 * Main Start
 *****************************************************************************/
int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s N (N must be a power of 2 and greater than 0.)\n", argv[0]);
        return 1;
    }

    FFT_N = atoi(argv[1]);
    if (FFT_N <= 0 || (FFT_N & (FFT_N - 1)) != 0) {
        fprintf(stderr, "N must be a power of 2 and greater than 0.\n");
        return 1;
    }

    prepare();
    time_gettime(&time_start);
    compute();
    time_gettime(&time_end);
    report();
    
    free(my);
    free(out);
    free(X_r);
    free(X_i);
    free(temp_r);
    free(temp_i);

    return 0;
}
/*****************************************************************************
 * End Main
 ***************************************************************************/

