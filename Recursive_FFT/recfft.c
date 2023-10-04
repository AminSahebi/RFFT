/*******************************************************************************
 *        Rec. FFT                                       			*
 *                                                                              *
 *  AUTHOR(s): Amin Sahebi							*
 *  VERSION: 0.0.1                                                              *
 *  DATE: 15/05/2019 (dd/mm/yy)                                                 *
 *  Based on the source code in https://www.math.wustl.edu                      *
 ******************************************************************************/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
/*******************************************************************************
 *Some Initial Definitions
 ******************************************************************************/
#define q       3	
#define N	(1<<q)	// N-point FFT 2^q points
double NUM_RET = 1000000;
typedef float real;
typedef struct{real Re; real Im;} complex;
complex v[N], scratch[N];

struct timeval tim;
double t1;
double t2;
#ifndef PI
# define PI	3.141592
#endif
//#define DEBUG

/*******************************************************************************
 *Recursive Fast Fourier Transform Body Function and Algorithm
 ******************************************************************************/
/*
 * This is the description of Recursive FFT algorithm 
 fft(v,N):
 [0] If N==1 then return.
 [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
 [2] Compute fft(ve, N/2);
 [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
 [4] Compute fft(vo, N/2);
 [5] For m = 0 to N/2-1, do [6] through [9]
 [6]   Let w.re = cos(2*PI*m/N)
 [7]   Let w.im = -sin(2*PI*m/N)
 [8]   Let v[m] = ve[m] + w*vo[m]
 [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */

void fft( complex *v, int n, complex *tmp )
{
	if(n<1) return;/* do nothing and return */

	if(n>1) {			
		int k,m;    complex z, w, *vo, *ve;
		ve = tmp; vo = tmp+n/2;
		for(k=0; k<n/2; k++) {
			ve[k] = v[2*k];
			vo[k] = v[2*k+1];
		}
		fft( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
		fft( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
		for(m=0; m<n/2; m++) {
			w.Re = cos(2*PI*m/(double)n);
			w.Im = -sin(2*PI*m/(double)n);
			z.Re = w.Re*vo[m].Re - w.Im*vo[m].Im;	/* Re(w*vo[m]) */
			z.Im = w.Re*vo[m].Im + w.Im*vo[m].Re;	/* Im(w*vo[m]) */
			v[m].Re = ve[m].Re + z.Re;
			v[m].Im = ve[m].Im + z.Im;
			v[m+n/2].Re = ve[m].Re - z.Re;
			v[m+n/2].Im = ve[m].Im - z.Im;
		}
	}
	return;
}


/*******************************************************************************
 *Prepare the input and other necessary elements
 ******************************************************************************/
void prepare()
{
	
#ifdef DEBUG
	printf("[INFO]: START prepare()\n");
#endif

	for (int i=0;i<N;i++){
		v[i].Re = 2*i+1; 
		v[i].Im = 0;
	}
#ifdef DEBUG
	for(int j=0;j<N;j++){
		printf("[INFO]: The inputs are %f , %f\n",v[j].Re,v[j].Im);
	}
#endif

}


/*******************************************************************************
 *Compute the FFT Algorithm, use iterations to normalize the execution time
 ******************************************************************************/

void compute()
{

	for (int k=0;k<NUM_RET;k++){
		fft( v, N, scratch );
	}
}


/*******************************************************************************
 *Report the Results
 ******************************************************************************/
void report(){

  printf("Execution Time for Recursive FFT %.6lf \n", (t2-t1)/NUM_RET);	
}

/*******************************************************************************
 * Main Body Function
 ******************************************************************************/
int main(void)
{
	t1 ,t2 =0.0;
	prepare();
	gettimeofday(&tim, NULL);
	t1=tim.tv_sec*1000000+(tim.tv_usec);
	compute();
 	gettimeofday(&tim, NULL);
	t2=tim.tv_sec*100000+(tim.tv_usec);
	report();

	return 0;
}
/*******************************************************************************
 *END.
 ******************************************************************************/
