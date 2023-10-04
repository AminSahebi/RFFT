#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//HOW to RUN
//gcc fftw.c -o fftw -lrt -lm -lfftw3


int main(int argc, char *argv[])
{
int NUM_ITERATION = 1000000;
int N = 8;
timeval tim;
double in[] = { 1, 3, 5, 7, 9, 11, 13, 15 };	/* Example input */

fftw_complex *out; /* Output */
fftw_plan p; /* Plan */

 /*
 * Size of output is (N / 2 + 1) because the other remaining items are
 * redundant, in the sense that they are complex conjugate of already
 * computed ones.
 *
 * CASE SIZE 6 (even):
 * [real 0][complex 1][complex 2][real 3][conjugate of 2][conjugate of 1]
 *
 * CASE SIZE 5 (odd):
 * [real 0][complex 1][complex 2][conjugate of 2][conjugate of 1]
 *
 * In both cases the items following the first N/2+1 are redundant.
 */
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

 /*
 * An fftw plan cares only about the size of in and out,
 * not about actual values. Can (and should) be re-used.
 */

gettimeofday(&tim,NULL);
double t1=tim.tv_sec*1000000+(tim.tv_usec);  
 
p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

gettimeofday(&tim, NULL);
double t2=tim.tv_sec*1000000+(tim.tv_usec);
printf("[INFO] FFTW Exectution Time %.6lf us\n", (t2-t1)/NUM_ITERATION);
 

  
 

//Execute the dft as indicated by the plan

 fftw_execute(p);

  
Print the N/2+1 complex values computed by the DFT function.

 int i;
for (i = 0; i < N ; i++) 
{
   printf("out[%d] = %f, %f\n", i, creal(out[i]), cimag(out[i]));
}

//Clean routine
  fftw_destroy_plan(p);
  fftw_free(out);

  return 1;
}
