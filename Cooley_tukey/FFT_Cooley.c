#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.141592655359
const float w = 0.7071067;

unsigned bit_reverse(register unsigned x){

	x=(((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555)<< 1));
	x=(((x & 0xcccccccc) >> 2) | ((x & 0x33333333)<< 2));
	x=(((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f)<< 4));
	x=(((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff)<< 8));

	return ((x >> 16) | (x << 16));
}

void fft_default(float real[], float imag[], unsigned n){

	//reordering
	unsigned log_n = __builtin_ctz(n);
	for(unsigned i=0;i<n;i++){
		unsigned j=bit_reverse(i)>>(sizeof(i)*8-log_n);
		if(i>j){

			float tmp_r[sizeof(i)*8-log_n];
			tmp_r[i] = real[i];
			real[i] = real[j];
			real[j] = tmp_r[i];

			float tmp_i[sizeof(i)*8-log_n];
			tmp_i[i] = imag[i];
			imag[i] = imag[j];
			imag[j] = tmp_i[i];

		}
	}

	for(unsigned i=1;i<n; i<<=1){
		for(unsigned k=0 ;k<i;k++){
			float w_real=cos(-2*PI*k/(i<<1));
			float w_imag=sin(-2*PI*k/(i<<1));

			for(unsigned j=0;j<n;j+=i<<1){

				float temp_real = real[j + k + i] * w_real - imag[j + k + i] * w_imag;
				float temp_imag = real[j + k + i] * w_imag + imag[j + k + i] * w_real;
				real[j + k + i] = real[j + k] - temp_real;
				imag[j + k + i] = imag[j + k] - temp_imag;
				real[j + k] += temp_real;
				imag[j + k] += temp_imag;

			}

		}

	}
}

	int main()
	{

		int FFT_N= 8;
		int number_of_samples = 1;
		float real[number_of_samples][FFT_N];
		float imag[number_of_samples][FFT_N];

		real[0][0] = 1; real[0][1]=3; real[0][2]= 5; real[0][3]=7; real[0][4]=9; real[0][5]=11; real[0][6]=13; real[0][7]=15;
		struct timeval tim;
		gettimeofday(&tim, NULL);
		double t3=tim.tv_sec+(tim.tv_usec);
		for (int k=0;k<100000;k++){
			for(int j=0; j<1;j++){
				fft_default(real[j],imag[j],8);
			}
		}
		gettimeofday(&tim, NULL);
		double t4=tim.tv_sec+(tim.tv_usec);

		printf("%.6lf Micro Seconds Elapsed for Default Algorithm without considering Data transfer\n", (t4-t3)/100000.0);
		return 0;
	}
