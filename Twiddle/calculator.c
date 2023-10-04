#include <stdio.h>
#include <math.h>


#define pi 3.14159265359

int main(int argc, char **argv){

int k;
double Re = 0;
double Im = 0;
sscanf(argv[1], "%d", &k);

for (int i=0; i<= k; i++){
Re = cos((2*pi*i)/k);
Im = -sin((2*pi*i)/k);
printf("the twiddle factor for  n=%d ,k=%d is %.12f  %.12fj \n",i,k, Re, Im);
}
return 0;

}
