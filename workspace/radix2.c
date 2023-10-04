/******************************************************************************
 *        Real-Valued FFT                                       	      *
 *                                                                            *
 *  AUTHOR(s): Amin Sahebi, Lorenzo Verdoscia, Roberto Giorgi - UNISI , CNR   *
 *  VERSION: 0.0.1                                                            *
 *  DATE: 01/05/2019 (dd/mm/yy)                                               *
 *                                                                            *
 *****************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

/*******************************************************************************
 *The Global Variables Definition
 ******************************************************************************/
#define PI 3.14159265359
#define MAXPOW 24
#define DEBUG
#define q 3
#define FFT_N (1<<q)
	int NUM_Iteration = 1;
	int    N;
	typedef struct {double R; double I;} complex;
	complex data[FFT_N];
	struct timeval tv;	
	double t1,t2;
	int pow_2[MAXPOW];
	int pow_4[MAXPOW];

void twiddle(complex *W, int N, double stuff)
{
	W->R=cos(stuff*2.0*PI/(double)N);
	W->I=-sin(stuff*2.0*PI/(double)N);
}


void bit_reverse_reorder(complex *W, int N)
{
	int bits, i, j, k;
	double tempr, tempi;

	for (i=0; i<MAXPOW; i++)
		if (pow_2[i]==N) bits=i;

	for (i=0; i<N; i++)
	{
		j=0;
		for (k=0; k<bits; k++)
			if (i&pow_2[k]) j+=pow_2[bits-k-1];

		if (j>i)  /** Only make "up" swaps */
		{
			tempr=W[i].R;
			tempi=W[i].I;
			W[i].R=W[j].R;
			W[i].I=W[j].I;
			W[j].R=tempr;
			W[j].I=tempi;
		}
	}
}

/** RADIX-2 FFT ALGORITHM */
void radix2(complex *data, int N)
{
	int    n2, k1, N1, N2;
	complex W, bfly[2];

	N1=2;
	N2=N/2;

	/** Do 2 Point DFT */
	for (n2=0; n2<N2; n2++)
	{
		/** Don't hurt the butterfly */
		twiddle(&W, N, (double)n2);
		bfly[0].R = (data[n2].R + data[N2 + n2].R);
		bfly[0].I = (data[n2].I + data[N2 + n2].I);
		bfly[1].R = (data[n2].R - data[N2 + n2].R) * W.R - 
			((data[n2].I - data[N2 + n2].I) * W.I); 
		bfly[1].I = (data[n2].I - data[N2 + n2].I) * W.R +
			((data[n2].R - data[N2 + n2].R) * W.I);

		/** In-place results */
		for (k1=0; k1<N1; k1++)
		{
			data[n2 + N2*k1].R = bfly[k1].R;
			data[n2 + N2*k1].I = bfly[k1].I;
		}
	}

	// Don't recurse if we're down to one butterfly//
	if (N2!=1)
		for (k1=0; k1<N1; k1++){
		*data = data[N2*k1];
		radix4(data, N2);
}
}


/*******************************************************************************
 *Prepare the input and other necessary elements
 ******************************************************************************/
void prepare(){

#ifdef DEBUG
	printf("[INFO]: START prepare()\n");
#endif
	N= FFT_N;
	
	for (int i=0;i<N;i++){
		data[i].R = 2.0*i+1.0; 
		data[i].I = 0.0;
	}

#ifdef DEBUG
	for(int j=0;j<FFT_N;j++){
		printf("[INFO]: The inputs are %f,%f\n",data[j].R,data[j].I);
	}
#endif
	/** Set up power of two arrays */
	pow_2[0]=1;
	for (int i=1; i<MAXPOW; i++)
		pow_2[i]=pow_2[i-1]*2;
	pow_4[0]=1;
	for (int i=1; i<MAXPOW; i++)
		pow_4[i]=pow_4[i-1]*4;

	}

/*******************************************************************************
 *Compute the FFT Algorithm, use iterations to normalize the execution time
 ******************************************************************************/
void compute()
{
	//for (int i=0 ; i< NUM_Iteration; i++){
	 radix2(data, FFT_N);
   // }
}

/*******************************************************************************
 *Report the Results
 ******************************************************************************/
void report(){
#ifdef DEBUG
	for (int i=0; i<N; i++)
		printf("The Results are %lf+j%lf\n",data[i].R,data[i].I);
#endif
	printf("Radix4 Execution time is %.6lf us\n", (t2-t1)/NUM_Iteration);
}

/*******************************************************************************
 * Main Body Function
 ******************************************************************************/

int main()
{
	
	prepare();
	
	gettimeofday(&tv,NULL);
	t1= (tv.tv_sec) + (tv.tv_usec);
	compute();
	
	gettimeofday(&tv,NULL);
	t2= (tv.tv_sec) + (tv.tv_usec);			
	bit_reverse_reorder(data, N);
	report();
return 0;

}



