#include <stdio.h>
#include <stdlib.h>
#include <time.h>



int main(){
struct timespec start,finish;

clock_gettime(CLOCK_REALTIME, &start); 

for(int i=0;i<100000000;i++){

	i++;
}

clock_gettime(CLOCK_REALTIME, &finish);

long long ns = finish.tv_nsec - start.tv_nsec;
printf("nanoseconds: %lld\n", ns);

return 0;
}


