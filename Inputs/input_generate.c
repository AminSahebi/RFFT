#include <stdio.h>
#include <stdlib.h>

#define NUM 4096
int main(void)
{

    char    str[NUM];
    FILE    *fptr;
    int     i;
    int     num;
    char    num2;
    i = 0;

    fptr = fopen("4096_input.txt", "w");
    if (fptr == NULL)
    {
        printf("ERROR Creating File!");
        exit(1);
    }
 	 for (int i = 0; i <NUM ; ++i){
  	fprintf(fptr, "%d", rand() % 10);
	fprintf(fptr, "\n");  /* Probably nice to make it a line. */
	}	
    	fclose(fptr);
    return 0;
}
