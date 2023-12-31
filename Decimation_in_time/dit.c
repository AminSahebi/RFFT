#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>


#define TRSIZ 8
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr

/******************************************************************************
 * FUNCTION: dittt()
 *
 * PURPOSE:
 * - Calculates Fast Fourier Transform of given data series using bit
 * reversal prior to FFT. This function is directly taken from the book
 * Numerical Recipes in C. The algorithm follows the butterfly diagram
 * 4 1 A BASIC DECIMATION-IN-TIME (DIT) ALGORITHM
 * of the radix-2 Decimation In Time (DIT) implementation of FFT, where
 * the bit reversal is mandatory. The function replaces the (input) array
 * data[1..2*nn] by its discrete Fourier transfor, if the parameter isign
 * is given as -1; or replaces data[1..2*nn] by nn times its inverse
 * discrete Fourier transform, if isign is given as 1. The array of data[]
 * is a complex array of length nn or, equivalently, a real array of
 * length 2*nn. nn MUST be an integer power of 2, but this is not checked!
 *
 * The needed increments of trigonometric functions are handled via
 * a clever identity using the fact that:
 *
 * exp(ia + ib) = exp(ia)*exp(ib) = exp(ia) + Z*exp(ia),
 *
 * where a is the starting angle and b is the increment. Solving for Z
 * gives:
 *
 * Z = exp(ib) - 1 = -(1 - cos(b)) + isin(b) = -2*sin^2(b/2) + isin(b).
 *
 * Then the increments can be calculated with respect to the zero-angle
 * set by wr (omega-real) and wi (omega-imaginary) by relation:
 *
 * (wr + wi) + Z*(wr + wi) = (wr + wi) + [-2*sin^2(b/2) + isin(b)](wr + wi).
 *
 * By setting wpr (omega-phase-real) = -2*sin^2(b/2) and wpi = sin(b),
 * one has
 *
 * (wr + wi) + (wpr + wpi)*(wr + wi) = exp(ia) + Z*exp(ia),
 *
 * where the real part is: wr + wpr*wr - wpi*wi
 * and the imaginary part is: wi + wpr*wi + wpi*wr.
 * The actual incremental parts here are:
 * wpr*wr - wpi*wi for real part and wpr*wi + wpi*wr for imaginary part.
 *
 *
 * INPUT:
 * - data[] = Array containing the data to be transformed
 * The transformed data is stored back to this array
 * so that the real and imaginary parts are following
 * each other --> size of array = 2*size of data
 * - nn = size of data = size of array/2, has to be a power of 2
 * - isign = if -1, calculates the FFT;
 * if 1, calculates the IFFT i.e. the inverse.
 *
 * OUTPUT: - (void)
 */

void dit()
{
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;
	int N = TRSIZ;
	int i = 0, j = 0, n = 0, k = 0, m = 0, isign = -1,istep,mmax;
	double data1[2*TRSIZ] = {1,0,3,0,5,0,7,0,9,0,11,0,13,0,15,0};
	double *data;
	data = &data1[0] - 1;
	n = N*2;
	j = 1;
	// do the bit-reversal
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j +=m;
	}
	// calculate the FFT
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
				j = i + mmax;
				tempr = wr*data[j] - wi*data[j+1];
				tempi = wr*data[j+1] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] = data[i] + tempr;
				data[i+1] = data[i+1] + tempi;
			}
			wtemp = wr;
			wr += wtemp*wpr - wi*wpi;
			wi += wtemp*wpi + wi*wpr;
		}
		mmax = istep;
	}

} 

int main(){

struct timeval tim;
double t1,t2;


	gettimeofday(&tim, NULL);
	t1=tim.tv_sec+(tim.tv_usec);
	for (int i=0;i<100000;i++){
		dit();
	}
	gettimeofday(&tim, NULL);
	t2=tim.tv_sec+(tim.tv_usec);
	// print the results
	printf("\nFourier components from the DIT algorithm:\n");
	printf("time Calculated for FFT is %.6f\n",(t2-t1)/100000.0);
	/*for (k = 0; k < 2*N; k +=2 )
		printf("\n%f %fi\n", data[k+1], data[k+2]);
*/
	return 0;
}
