#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
//#include <sys/time.h>
#include <time.h>
int i;
int N = 8;
//timeval tim;
double number_of_iteration = 1000;
double in[] = { 1, 3, 5, 7, 9, 11, 13, 15 };	/* Example input */
struct timespec start,finish;

int  main()
{

	fftw_plan p; /* Plan */
	fftw_complex *out; /* Output */
	clock_gettime(CLOCK_REALTIME, &start);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
	for(i=0; i<number_of_iteration ;i++){

		/*
		 *    * Size of output is (N / 2 + 1) because the other remaining items are
		 *    * redundant, in the sense that they are complex conjugate of already
		 *    * computed ones.
		 *    *
		 *    * CASE SIZE 6 (even):
		 *    * [real 0][complex 1][complex 2][real 3][conjugate of 2][conjugate of 1]
		 *    *
		 *    * CASE SIZE 5 (odd):
		 *    * [real 0][complex 1][complex 2][conjugate of 2][conjugate of 1]
		 *    *
		 *    * In both cases the items following the first N/2+1 are redundant. 
		 *    * An fftw plan cares only about the size of in and out,
		 *       * not about actual values. Can (and should) be re-used.
		 *          */

		p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
		/********************************************
		 *  Execute the dft as indicated by the plan
		 **********************************************/
		fftw_execute(p);

	}
	clock_gettime(CLOCK_REALTIME, &finish);

	fftw_destroy_plan(p);
	fftw_free(out);
	double  ns = finish.tv_nsec - start.tv_nsec;

	printf("time Elapsed by nanosecond %5.5f\n",ns/number_of_iteration);
	/************************************************************
	 ** Print the N/2+1 complex values computed by the DFT function.
	 ************************************************************/

	return 0;
}


