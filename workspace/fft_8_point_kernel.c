/*******************************************************************************
 *        Real-Valued FFT kernel based on 8_point                                      		*
 *                                                                              *
 *  AUTHOR(s): Amin Sahebi, Lorenzo Verdoscia, Roberto Giorgi - UNISI , CNR     *
 *  VERSION: 0.0.1                                                              *
 *  DATE: 10/09/2019 (dd/mm/yy)                                                 *
 *                                                                              *
 ******************************************************************************/
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

/*******************************************************************************
 *The Global Variables Definition
 ******************************************************************************/
typedef unsigned char byte;
const float w = 0.70710;
double t1,t2;
float X_i[8] , temp_i[8];	
float *out_r,*out_i,*data,*my,*X_r,*temp_r;
struct timeval tim;
int Num_Iterations = 1;
int FFT_N = 16;

//#define DEBUG 

/*******************************************************************************
 *bit-reversal function to find even and odd numbers
 ******************************************************************************/


//inline void bit_reverse(int *ar, int size){
/*static const unsigned char BitReverseTable256[] =// {0x00, 0x08, 0x04, 0x0C, 0x02, 0x0A, 0x06, 0x0E, 0x01, 0x09, 0x05, 0x0D, 0x03, 0x0B, 0x07, 0x0F };
  {
  0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2, 
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
  };

  int bit_reverse(int x){

  return ((BitReverseTable256[x & 0xff] << 24) | 
  (BitReverseTable256[(x >> 8) & 0xff] << 16) | 
  (BitReverseTable256[(x >> 16) & 0xff] << 8) |
  (BitReverseTable256[(x >> 24) & 0xff]));
  }
  */



/*static inline __attribute__ ((always_inline))
  unsigned int bit_reverse (unsigned int x)
  {
  register unsigned int y = 0x55555555;
  x = (((x >> 1) & y) | ((x & y) << 1));
  y = 0x33333333;
  x = (((x >> 2) & y) | ((x & y) << 2));
  y = 0x0f0f0f0f;
  x = (((x >> 4) & y) | ((x & y) << 4));
  y = 0x00ff00ff;
  x = (((x >> 8) & y) | ((x & y) << 8));
  return((x >> 16) | (x << 16));
  }
  */
unsigned bit_reverse(register int x){ //low-Memory

	x=(x & 0xaaaaaaaa) >> 1 | (x & 0x55555555) << 1;
	x=(x & 0xcccccccc) >> 2 | (x & 0x33333333) << 2;
	x=(x & 0xf0f0f0f0) >> 4 | (x & 0x0f0f0f0f) << 4; 
	x=(x & 0xff00ff00) >> 8 | (x & 0x00ff00ff) << 8;

	return ((x >> 16 ) | (x << 16 ));
}
/*
   inlie byte bit_reverse(byte c); __attribute__ ((always_inline))
   inline byte bit_reverse(byte c)
   {
   c = ((c>>1)&0x55)|((c<<1)&0xAA);
   c = ((c>>2)&0x33)|((c<<2)&0xCC);
   return ((c>>4) | (c<<4));
   }
   */
/*******************************************************************************
 *Real Valued Fast Fourier Transform Body Function
 ******************************************************************************/

void fft_core(float my[],int n)
{ 
	//struct C_Num X[8], tmp[8];
	float a_r, a_i;
	unsigned log_n = __builtin_ctz(n);

	for(int i=0;i<n;i++)
	{
		unsigned j=bit_reverse(i)>>(sizeof(i)*8-log_n); 	  
		if(i>j){
			float tmp_r[sizeof(i)*8-log_n];  //re-ordering
			tmp_r[i] = my[i];
			my[i] = my[j];
			my[j] = tmp_r[i];
		}
	}

#ifdef DEBUG
	for(int i=0;i<n;i++)
	{
		printf("[INFO]: The re-ordered elements are...%f\n",my[i]);
	}
#endif
	memset(X_i, 0, 8 * sizeof(float));
	memset(temp_i, 0, 10 * sizeof(float));

	memcpy(temp_r,my, 8 * sizeof(float));
	memset(temp_i, 0, 10*sizeof(float));
#ifdef DEBUG
	for(int i=0;i<n;i++)
	{
		printf("[INFO]: The TEMP elements are...%f,%f\n",temp_r[i],temp_i[i]);
	}
#endif

	/*for (int n=0; n<8; n++)//bottleneck changed with MEMCPY and MEMSET
	  {
	  X[n].i =0.0;
	  tmp[n].i =0.0;
	  tmp[n].R = my[n];
	  }*/

	for(int n=0; n<4; n++)
	{
		X_r[2*n] = temp_r[2*n] + temp_r[2*n+1];
		X_r[2*n+1] = temp_r[2*n] - temp_r[2*n+1];
	}

	memcpy(temp_r, X_r, 8 * sizeof(float));
	memcpy(temp_i ,X_i, 8 * sizeof(float));
	/*for(int n=0; n<8; n++) //bottleneck
	  {
	  tmp[n].R = X[n].R;
	  tmp[n].i = X[n].i;
	  }*/

	temp_i[3] = - temp_r[3]; temp_r[3] = 0;
	temp_i[7] = - temp_r[7]; temp_r[7] = 0;

	for(int n=0; n<2; n++)
	{
		X_r[4*n] = temp_r[4*n] + temp_r[4*n+2];
		X_r[4*n+2] = temp_r[4*n] - temp_r[4*n+2];
		X_r[4*n+1] = temp_r[4*n+1] + temp_r[4*n+3];
		X_i[4*n+1] = temp_i[4*n+1] + temp_i[4*n+3];
		X_r[4*n+3] = temp_r[4*n+1] - temp_r[4*n+3];
		X_i[4*n+3] = temp_i[4*n+1] - temp_i[4*n+3];
	}


	memcpy(temp_r, X_r, 8 * sizeof(float));
	memcpy(temp_i, X_i, 8 * sizeof(float));
	/*for(int n=0; n<8; n++) //bottleneck
	  {
	  tmp[n].R = X[n].R;
	  tmp[n].i = X[n].i;
	  } */


	temp_i[6] = -temp_r[6]; temp_r[6] = 0;
	a_r = w*(temp_r[5] - temp_i[7]);
	a_i = -w*(temp_r[5] + temp_i[7]);
	temp_r[5]=a_r; temp_i[5]=a_i;
	temp_r[7]=a_r; temp_i[7]=a_i;

	for(int n=0; n<4; n++)
	{
		X_r[n] = temp_r[n] + temp_r[n+4];
		X_i[n] = temp_i[n] + temp_i[n+4];
		X_r[n+4] = temp_r[n] - temp_r[n+4];
		X_i[n+4] = temp_i[n] - temp_i[n+4];
	}
#ifdef DEBUG

	for(int i=0;i<n;i++)
	{
		printf("[INFO]: RFFT calculated outputs are...%f+j%f\n",X_r[i], X_i[i]);
	}
#endif


}//end of FFT Function

/*******************************************************************************
 *Prepare the input and other necessary elements
 ******************************************************************************/
void prepare()
{
#ifdef DEBUG
	printf("[INFO]: START prepare()\n");
#endif
	data = (float*) malloc(FFT_N*sizeof(float));
	out_r =(float*) malloc(FFT_N*sizeof(float));
	out_i =(float*) malloc(FFT_N*sizeof(float));
	my = (float*) malloc(8*sizeof(float));
	for (int i = 0; i < FFT_N; i++) 
	{
		data[i] = 2.0*i+ 1;	
	}

	X_r = (float*)malloc(sizeof(float));
	temp_r = (float*)malloc(sizeof(float));
#ifdef DEBUG
	for(int j=0;j<FFT_N;j++)
	{
		printf("[INFO]: The inputs in the memory are %f\n",data[j]);
	}
#endif

}

/*******************************************************************************
 *Compute the FFT Algorithm, use iterations to normalize the execution time
 ******************************************************************************/
int compute()
{
		for (int i = 0; i < Num_Iterations; i++)
		{ 
			for (int k = 0 ; k < FFT_N ; k+=8)
			{
					for (int j =0 ; j < 8 ; j++)
					{
					my[j] = data[j+k];
					}
						fft_core(my,8);
							for (int l=0 ; l< FFT_N; l++){
								out_r[l+k]=X_r[l];
								out_i[l+k]=X_i[l];
						}
					}	
				}
	return 0;
}


		/*******************************************************************************
		 *Report the Results
		 ******************************************************************************/

		void report(){

			printf("%.6lf [info]Execution Time Proposed Algorithm \n", (t2-t1)/Num_Iterations);

		}
		/*****************************************************************************
		 *  Main Start
		 *****************************************************************************/
		int main(int argc, char *argv[])
		{

			prepare();
			gettimeofday(&tim, NULL);
			t1=tim.tv_sec*1000000+(tim.tv_usec);

			compute();

			gettimeofday(&tim, NULL);
			t2=tim.tv_sec*1000000+(tim.tv_usec);

			report();
			return 0;
		}
		/*****************************************************************************
		 * End Main
		 ***************************************************************************/
