#include <stdio.h>
#include "math.h"
#include <stdlib.h>

// lokalnie
#include "sources/nr.h"
#include "sources/nrutil.h"
#include "sources/nrutil.c"
#include "sources/gauleg.c"
#include "sources/gaulag.c"
#include "sources/gauher.c"
#include "sources/gammln.c"

/*
// z taurusa
#include "/opt/NR/numerical_recipes.c/nr.h"
#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gauleg.c"
#include "/opt/NR/numerical_recipes.c/gaulag.c"
#include "/opt/NR/numerical_recipes.c/gauher.c"
#include "/opt/NR/numerical_recipes.c/gammln.c"
*/

////////////////////////////////////////////////////////////////////// AD 1 - g(x_i)
float g1(float x)
{
	return 1/(x*sqrt(pow(x,2)-1));
}

////////////////////////////////////////////////////////////////////// AD 2a - g(x_i)
float g2a(float x)
{
	return log(fabs(x));
}

////////////////////////////////////////////////////////////////////// AD 2b - g(x_i)
float g2b(float x)
{
	return log(x)*exp(-pow(x,2));
}

////////////////////////////////////////////////////////////////////// AD 3 - g(x_i)
float g3(float x)
{
	return sin(2*x)*exp(-2*x);
}

/*
	FILE *file1;
	file1 = fopen("file1.txt","w");
	fclose(file1);
*/

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
	float x_min, x_max, C_dokl, I;
	int n, n_max = 100;
	float *x = (float*)malloc((n_max+1)*sizeof(float));
	float *w = (float*)malloc((n_max+1)*sizeof(float));
////////////////////////////// AD 1 - GAUSS-LEGANDRE
	FILE *gauleg1;
	gauleg1 = fopen("gauleg1.txt","w");
	x_min = 1;
	x_max = 2;
	C_dokl = M_PI/3;
	for(n=2; n<=100; ++n)
	{
		gauleg(x_min,x_max,x,w,n);
	
		//for(int i=1; i<=n; ++i) printf("w_%d = %lf\tx_%d = %lf\n", i, w[i], i, x[i]);
	
		I = 0;
		for(int i=1; i<=n; ++i) I += w[i]*g1(x[i]);
	
		//printf("I = %f\n", I);
		fprintf(gauleg1, "%d\t%f\n",n,fabs(C_dokl - I));
	}
	fclose(gauleg1);
////////////////////////////// AD 2a - GAUSS-HERMITE
	FILE *gauher2a;
	gauher2a = fopen("gauher2a.txt","w");
	x_min = 0;
	x_max = INFINITY;
	C_dokl = -0.8700577;
	for(n=2; n<=100; n+=2)
	{
		gauher(x,w,n);
	
		//for(int i=1; i<=n; ++i) printf("w_%d = %lf\tx_%d = %lf\n", i, w[i], i, x[i]);
	
		I = 0;
		for(int i=1; i<=n; ++i) I += w[i]*g2a(x[i]);
		I = I/2;
	
		//printf("I = %f\n", I);
		fprintf(gauher2a,"%d\t%f\n",n,fabs(C_dokl - I));
	}
	fclose(gauher2a);
////////////////////////////// AD 2b - GAUSS-LEGENDRE
	FILE *gauleg2b;
	gauleg2b = fopen("gauleg2b.txt","w");
	x_min = 0;
	x_max = 5;
	C_dokl = -0.8700577;
	for(n=2; n<=100; ++n)
	{
		gauleg(x_min,x_max,x,w,n);
	
		//for(int i=1; i<=n; ++i) printf("w_%d = %lf\tx_%d = %lf\n", i, w[i], i, x[i]);
	
		I = 0;
		for(int i=1; i<=n; ++i) I += w[i]*g2b(x[i]);
	
		//printf("I = %f\n", I);
		fprintf(gauleg2b, "%d\t%f\n",n,fabs(C_dokl - I));
	}
	fclose(gauleg2b);
////////////////////////////// AD 3 - GAUSS-LAGUERRE
	FILE *gaulag3;
	gaulag3 = fopen("gaulag3.txt","w");
	x_min = 0;
	x_max = INFINITY;
	C_dokl = 2./13;
	for(n=2; n<=10; ++n)
	{
		gaulag(x,w,n,0);
		//for(int i=1; i<=n; ++i) printf("w_%d = %lf\tx_%d = %lf\n", i, w[i], i, x[i]);
	
		I = 0;
		for(int i=1; i<=n; ++i) I += w[i]*g3(x[i]);
	
		//printf("I = %f\n", I);
		fprintf(gaulag3,"%d\t%f\n",n,fabs(C_dokl - I));
	}
	fclose(gaulag3);
	free(x);
	free(w);
	return 0;
}

