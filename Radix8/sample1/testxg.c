/* test of fft*g.c */

#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#define MAX(x,y) ((x) > (y) ? (x) : (y))

/* random number generator, 0 <= RND < 1 */
#define RND(p) ((*(p) = (*(p) * 7141 + 54773) % 259200) * (1.0 / 259200.0))

#ifndef NMAX
#define NMAX 8192
#define NMAXSQRT 64
#endif

void cdft(int, int, double *, int *, double *);
void rdft(int, int, double *, int *, double *);
void ddct(int, int, double *, int *, double *);
void ddst(int, int, double *, int *, double *);
void dfct(int, double *, double *, int *, double *);
void dfst(int, double *, double *, int *, double *);
void putdata(int nini, int nend, double *a);
double errorcheck(int nini, int nend, double scale, double *a);


int main()
{	int num_iteration = 10000;
	int n, ip[NMAXSQRT + 2];
	double a[NMAX + 1], w[NMAX * 5 / 4], t[NMAX / 2 + 1], err;
	struct timeval tim;
	printf("data length n=? (must be 2^m)\n");
	scanf("%d", &n);
	ip[0] = 0;

	/* check of CDFT */
	/*    putdata(0, n - 1, a);
	      cdft(n, 1, a, ip, w);
	      cdft(n, -1, a, ip, w);
	      err = errorcheck(0, n - 1, 2.0 / n, a);
	      printf("cdft err= %g \n", err);
	      */
	/* check of RDFT */
	putdata(0, n - 1, a);
	gettimeofday(&tim , NULL);
	double t1=(tim.tv_sec) * 1000000 + (tim.tv_usec) ;
	for(int i=0; i< num_iteration; i++){
		rdft(n, 1, a, ip, w);
	}
	gettimeofday(&tim , NULL);
	double t2=(tim.tv_sec) * 1000000 + (tim.tv_usec) ;
	for(int i=0; i< num_iteration; i++){
		rdft(n, -1, a, ip, w);
	}
	err = errorcheck(0, n - 1, 2.0 / n, a);
	printf("rdft err= %g \n", err);
	printf("Execution time is %f\n",(t2-t1)/num_iteration);

	/* check of DDCT */
	/*    putdata(0, n - 1, a);
	      ddct(n, 1, a, ip, w);
	      ddct(n, -1, a, ip, w);
	      a[0] *= 0.5;
	      err = errorcheck(0, n - 1, 2.0 / n, a);
	      printf("ddct err= %g \n", err);
	      */
	/* check of DDST */
	/*    putdata(0, n - 1, a);
	      ddst(n, 1, a, ip, w);
	      ddst(n, -1, a, ip, w);
	      a[0] *= 0.5;
	      err = errorcheck(0, n - 1, 2.0 / n, a);
	      printf("ddst err= %g \n", err);
	      */
	/* check of DFCT */
	/*    putdata(0, n, a);
	      a[0] *= 0.5;
	      a[n] *= 0.5;
	      dfct(n, a, t, ip, w);
	      a[0] *= 0.5;
	      a[n] *= 0.5;
	      dfct(n, a, t, ip, w);
	      err = errorcheck(0, n, 2.0 / n, a);
	      printf("dfct err= %g \n", err);
	      */
	/* check of DFST */
	/*    putdata(1, n - 1, a);
	      dfst(n, a, t, ip, w);
	      dfst(n, a, t, ip, w);
	      err = errorcheck(1, n - 1, 2.0 / n, a);
	      printf("dfst err= %g \n", err);
	      */
	return 0;
}


void putdata(int nini, int nend, double *a)
{
	int j, seed = 0;

	for (j = nini; j <= nend; j++) {
		a[j] = RND(&seed);
	}
}


double errorcheck(int nini, int nend, double scale, double *a)
{
	int j, seed = 0;
	double err = 0, e;

	for (j = nini; j <= nend; j++) {
		e = RND(&seed) - a[j] * scale;
		err = MAX(err, fabs(e));
	}
	return err;
}

