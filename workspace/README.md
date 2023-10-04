-- RFFT The C files here, prepared and used in this article, to use our proposed algorithm please cite the paper below,

L. Verdoscia, A. Sahebi and R. Giorgi, "A Data-Flow Methodology for Accelerating FFT," 2019 8th Mediterranean Conference on Embedded Computing (MECO), Budva, Montenegro, 2019, pp. 1-4. doi: 10.1109/MECO.2019.8760044

This is the C/C++ code files to show our proposed new FFT algorithm implementation to test and verify the methodology latency. moreover, we include some codes from FFTW Ver. 3.3.8, fftw_plan_dft_r2c_1d plan, and the default implementation based on Cooley-Tukey FFT.

To build the FFTW code you have to do following command, "g++ fft_real_valued.cpp -o fft_real_value -lfftw3 -lm" for any furthur description about how to setup FFTW in your machine please visit the www.fftw.org and you can read more about at: http://www.fftw.org/fftw3.pdf
