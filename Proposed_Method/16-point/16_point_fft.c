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
struct timeval tim;
int Num_Iterations = 100000;
int FFT_N = 16;
#define pi 3.14159265359
//#define DEBUG 

/*******************************************************************************
 *bit-reversal function to find even and odd numbers
 ******************************************************************************/
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
		printf("[INFO]: The re-ordered elements are...%f\n",my[i]);
#endif

	for (int k = 0 ; k < FFT_N ; k+=8)//begining of 8-RFFT
	{
		memset(temp_r, 0, 8 * sizeof(DATA));
		memset(X_i, 0, 8 * sizeof(DATA));
		memset(temp_i, 0, 8 * sizeof(DATA));
		memset(X_r, 0, 8 * sizeof(DATA));
		for (int j =0 ; j < 8 ; j++){
			temp_r[j] = my[j+k];
		}
/*#ifdef DEBUG 
	for(int i =0; i<8;i++)
	printf("here temp_r is %f\n",temp_r[i]);
#endif*/
		for(int n=0; n<4; n++)
		{
			X_r[2*n] = temp_r[2*n] + temp_r[2*n+1];
			X_r[2*n+1] = temp_r[2*n] - temp_r[2*n+1];
		}

		memcpy(temp_r, X_r, 8 * sizeof(DATA));
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

		memcpy(temp_r, X_r, 8 * sizeof(DATA));
		memcpy(temp_i, X_i, 8 * sizeof(DATA));

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
		printf("[INFO]: RFFT calculated outputs out_r , out_i are...%f+j%f\n",out_r[i], out_i[i]);
#endif
        for(int i=8; i< FFT_N ; i++){//Multiplication
                out_r[i]=(out_r[i] * W_r[i-8])-(out_i[i] * W_i[i-8]);
                out_i[i]=(out_r[i] * W_i[i-8])+(out_i[i] * W_r[i-8]);
        }
#ifdef DEBUG
                for(int i=8; i< FFT_N;i++)
                printf("the outputs after MULT are X[%d]  %f %f \n",i, out_r[i] , out_i[i]);
#endif
                for(int i=0; i< FFT_N/2;i++){
                R_r[i]=out_r[i]+out_r[i+8];
                R_i[i]=out_i[i]+out_i[i+8];
                R_r[i+8]=-out_r[i+8]+out_r[i];
                R_i[i+8]=-out_i[i+8]+out_i[i];
                }
	
                for(int i=0; i< FFT_N;i++){
		R_i[i]+=(R_r[i]*W_i[i]);
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
	out_r =(DATA*) calloc(FFT_N,sizeof(DATA));
	out_i =(DATA*) calloc(FFT_N,sizeof(DATA));
	my = (DATA*) calloc(FFT_N,sizeof(DATA));
	data = (DATA*) calloc(FFT_N,sizeof(DATA));
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
	for (int i=0; i<FFT_N; i++){
	    W_r[i] = cos((2*pi*i)/(DATA)16);
	}
	for(int i=0; i<FFT_N;i++){
	    W_i[i] = -sin((2*pi*i)/(DATA)16);
	}
#ifdef DEBUG
	for(int j=0;j<FFT_N;j++)
		printf("[INFO]: The inputs in the memory are %f\n",my[j]);
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
/* 	double sum_input = 0; 
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
        printf( Err < 0.0001 ? "*** SUCCESS ***\n" : "*** FAIL ***\n");*/
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
