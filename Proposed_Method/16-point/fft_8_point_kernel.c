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
#define DATA float
typedef unsigned char byte;
float w = 0.707106781187;
double t1=0;
double t2=0;

float *out_r,*out_i,*data,*my,*X_r,*temp_r,*X_i, *temp_i, *R_r, *R_i, *W_r, *W_i;
float *Even_r, *Even_i, *Odd_r, *Odd_i;
float *MULT_r,*MULT_i;
struct timeval tim;
int Num_Iterations = 1;
int FFT_N = 16;
#define pi 3.14159265359
#define DEBUG 

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
/*******************************************************************************
 *Real Valued Fast Fourier Transform Body Function
 ******************************************************************************/
void fft_core(float my[],int n)
{ 
	float a_r, a_i;
	unsigned log_n = __builtin_ctz(FFT_N);

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

	for (int k = 0 ; k < FFT_N ; k+=8)//begining of 8-RFFT
	{

		memset(temp_r, 0, 8 * sizeof(float));
		memset(X_i, 0, 8 * sizeof(float));
		memset(temp_i, 0, 8 * sizeof(float));
		memset(X_r, 0, 8 * sizeof(float));
		for (int j =0 ; j < 8 ; j++){
			temp_r[j] = my[j+k];
		}
/*#ifdef DEBUG 
	printf("-------------------------------------\n");
	for(int i =0; i<8;i++)
	printf("here temp_r is %f\n",temp_r[i]);
#endif*/
	
		for(int n=0; n<4; n++)
		{
			X_r[2*n] = temp_r[2*n] + temp_r[2*n+1];
			X_r[2*n+1] = temp_r[2*n] - temp_r[2*n+1];
		}

		memcpy(temp_r, X_r, 8 * sizeof(float));
		//memcpy(temp_i ,X_i, 8 * sizeof(float));
		
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
		for (int l=0 ; l< FFT_N/2; l++){
			out_r[l+k]=X_r[l];
			out_i[l+k]=X_i[l];
		}

	}//end of 8-RFFT
#ifdef DEBUG
		for(int i=0;i<FFT_N;i++)
		{
		printf("[INFO]: RFFT calculated outputs out_r , out_i are...%f+j%f\n",out_r[i], out_i[i]);
		}
#endif
/*
	for(int k=0; k<FFT_N/2;k++){
		Even_r[k]=out_r[2*k];
		Even_i[k]=out_i[2*k];
		Odd_r[k]=out_r[2*k+1];
		Odd_i[k]=out_i[2*k+1];
	}

#ifdef DEBUG
		for(int i=0;i<FFT_N/2;i++)
		{
		printf("%f+j%f  -------    %f+j%f\n",Even_r[i], Even_i[i], Odd_r[i], Odd_i[i]);
		}
#endif
	
	for(int i=0; i< FFT_N/2 ; i++){//Multiplication
		MULT_r[i]=(Odd_r[i] * W_r[i])-(Odd_i[i] * W_i[i]);
		MULT_i[i]=(Odd_r[i] * W_i[i])+(Odd_i[i] * W_r[i]);
	}	
#ifdef DEBUG
		for(int i=0; i< FFT_N/2;i++)	
		printf("the outputs after MULT are X[%d] %f  %f\n",i, MULT_r[i] , MULT_i[i]);
#endif

	for(int i=0; i< FFT_N/2 ; i++){//finalize
		R_r[i]=MULT_r[i]+Even_r[i];
		R_i[i]=MULT_i[i]+Even_i[i];
		R_r[i+FFT_N/2]=Even_r[i]-MULT_r[i];
		R_i[i+FFT_N/2]=Even_i[i]-MULT_i[i];
		}
#ifdef DEBUG
		for(int i=0; i<FFT_N;i++)
	printf("the 16 - point results %f  %f\n",R_r[i],R_i[i]);
#endif
*/

        for(int i=0; i< FFT_N/2 ; i++){//Multiplication
                MULT_r[i]=(out_r[2*i+1] * W_r[i])-(out_i[2*i+1] * W_i[i]);
                MULT_i[i]=(out_r[2*i+1] * W_i[i])+(out_r[2*i+1] * W_r[i]);
        }
#ifdef DEBUG
                for(int i=0; i< FFT_N/2;i++)
                printf("the outputs after MULT are X[%d]  %f,%f \n",i, MULT_r[i] , MULT_i[i]);
#endif
                for(int i=0; i< FFT_N/2;i++){
                R_r[i]=out_r[i]+out_r[i+8];
                R_i[i]=out_i[i]+out_i[i+8];
                }
                for(int i=0; i<FFT_N/2;i++)
                {
                R_r[i+8]=-out_r[i+8]+out_r[i];
                R_i[i+8]=-out_i[i+8]+out_i[i];
                }
#ifdef DEBUG
                for(int i=0; i<FFT_N;i++)
        printf("the 16 - point results %f  %f\n",R_r[i],R_i[i]);
#endif


//
}//end of FFT Function

/*******************************************************************************
 *Prepare the input and other necessary elements
 ******************************************************************************/
void prepare()
{
#ifdef DEBUG
	printf("[INFO]: START prepare()\n");
#endif
	//data = (float*) calloc(FFT_N,sizeof(float));
	R_r= (DATA*) calloc(FFT_N, sizeof(DATA));
	R_i= (DATA*) calloc(FFT_N, sizeof(DATA));
	MULT_r= (DATA*) calloc(FFT_N/2, sizeof(DATA));
	MULT_i= (DATA*) calloc(FFT_N/2, sizeof(DATA));
	Even_i= (DATA*) calloc(FFT_N/2, sizeof(DATA));
	Even_r= (DATA*) calloc(FFT_N/2, sizeof(DATA));
	Odd_r= (DATA*) calloc(FFT_N/2, sizeof(DATA));
	Odd_i= (DATA*) calloc(FFT_N/2, sizeof(DATA));
	out_r =(DATA*) calloc(FFT_N,sizeof(DATA));
	out_i =(DATA*) calloc(FFT_N,sizeof(DATA));
	my = (DATA*) calloc(FFT_N,sizeof(DATA));
	X_r = (DATA*)calloc(8,sizeof(DATA));
	temp_r = (DATA*)calloc(8,sizeof(DATA));
	X_i = (DATA*)calloc(8,sizeof(DATA));
	temp_i = (DATA*)calloc(8,sizeof(DATA));
	W_i = (DATA*)calloc(8,sizeof(DATA));
	W_r=(DATA*)calloc(8,sizeof(DATA));
	for (int i = 0; i < FFT_N; i++) 
	{
		my[i] = (DATA)(2.0*i+ 1);	
	}
	for (int i=0; i<8; i++){
	    W_r[i] = cos((2*pi*i)/16);
	}
	for(int i=0; i<8;i++){
	    W_i[i] = -sin((2*pi*i)/16);
	}
#ifdef DEBUG
	for(int j=0;j<FFT_N;j++)
	{
		printf("[INFO]: The inputs in the memory are %f\n",my[j]);
	}
#endif

}

/*******************************************************************************
 *Compute the FFT Algorithm, use iterations to normalize the execution time
 ******************************************************************************/
void compute()
{
	for (int i = 0; i < Num_Iterations; i++)
	{ 
		fft_core(my,FFT_N);
	}
}


/*******************************************************************************
 *Report the Results
 ******************************************************************************/

void report(){

 	double sum_input = 0; 
        double sum_output = 0;

        for (int i=0;i<FFT_N;i++)
        {
                sum_input+=my[i]*my[i];
                sum_output+=(R_r[i]*R_r[i] + R_i[i]*R_i[i]);
        }
        double Err = sum_input-(sum_output/FFT_N);
        printf("the abs of inputs are %f\n",sum_input);
        printf("the abs of out puts are %f\n",sum_output/FFT_N);
        printf("The Err is %.6lf\n",sum_input-(sum_output/FFT_N));
        printf( Err < 0.0001 ? "*** SUCCESS ***\n" : "*** FAIL ***\n");
        printf("[info]Execution Time %.6lf\n", (t2-t1)/Num_Iterations);
	

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
