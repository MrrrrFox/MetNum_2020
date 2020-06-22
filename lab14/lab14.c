#include <stdio.h>
#include "math.h"
#include <stdlib.h>

////////////////////////////////////////////////////////////////////// wypisanie w kolumnach
void print(double tab[], int rozm, int start)
{
	for(int i=start; i<rozm-3; ++i) printf("%lf\t%lf\t%lf\t%lf\n", tab[i], tab[i+1], tab[i+2], tab[i+3]);
}

////////////////////////////////////////////////////////////////////// wypisanie w kolumnach do pliku
void fprint(FILE * file, double tab[], int rozm, int start)
{
	for(int i=start; i<rozm-3; ++i) fprintf(file, "%lf\t%lf\t%lf\t%lf\n", tab[i], tab[i+1], tab[i+2], tab[i+3]);
}

////////////////////////////////////////////////////////////////////// skalar
double skalar(double tab[], double tab2[])
{
	double sum = 0;
	for(int i=0; i<sizeof(tab)/sizeof(tab[0]); ++i) sum +=tab[i]*tab2[i];
	return sum;
}

////////////////////////////////////////////////////////////////////// normal euklidesowa
double norma(double tab[])
{
	return sqrt(skalar(tab,tab));
}

////////////////////////////////////////////////////////////////////// gen 1
double gen_1()
{
	static long int x=10;
	int a=17;
	int c=0;
	long int m=pow(2,13)-1;
	x = (a*x+c)%m;
	return x/(m+1.0);
}

////////////////////////////////////////////////////////////////////// gen 2
double gen_2()
{
	static long int x=10;
	int a=85;
	int c=0;
	long int m=pow(2,13)-1;
	x = (a*x+c)%m;
	return x/(m+1.0);
}

////////////////////////////////////////////////////////////////////// gen 3
double gen_3()
{
	static long long int x1=10,x2=10,x3=10;
	long long int m = (long long)pow(2,32) - 5;
	long long int x3_tmp = x3;
	x3 = (1176*x3 + 1476*x2 + 1776*x1 )%m;
	long long int x2_tmp = x2;
	x2 = x3_tmp;
	x1 = x2_tmp;
	return x3/(m+1.0);
}

/*
	FILE *file1;
	file1 = fopen("file1.txt","w");
	fclose(file1);
*/

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
// pliki
	FILE *gen1, *gen2, *gen3, *sfera, *kula, *gest1, *gest2, *gest3;
	gen1 = fopen("gen1.txt","w");
	gen2 = fopen("gen2.txt","w");
	gen3 = fopen("gen3.txt","w");
	sfera = fopen("sfera.txt","w");
	kula = fopen("kula.txt","w");
	gest1 = fopen("gest1.txt","w");
	gest2 = fopen("gest2.txt","w");
	gest3 = fopen("gest3.txt","w");
	FILE **gest = {gest1,gest2,gest3};
	int N_tab[] = {2000,pow(10,4),pow(10,7)};
	
	double *x1, *x2, *x3;
	double **r, *x, *y, *z;
	double **R;
	
	for(int licznik=0;licznik<sizeof(N_tab)/sizeof(N_tab[0]); ++licznik)
	{
		int N = N_tab[licznik];
	// GENERATOR 1.1 + 1.2
		x1 = (double*)malloc((N+1)*sizeof(double));
		x1[0] = 10;
		x2 = (double*)malloc((N+1)*sizeof(double));
		x2[0] = 10;
		for(int i=1; i<=N; ++i)
		{
			x1[i] = gen_1();
			x2[i] = gen_2();
		}
		if(licznik == 0) fprint(gen1,x1,N+1,1);
		if(licznik == 0) fprint(gen2,x2,N+1,1);
	// GENERATOR 1.3
		x3 = (double*)malloc((N+3)*sizeof(double));
		x3[0] = x3[1] = x3[2] = 10;
		for(int i=3; i<N+3; ++i) x3[i] = gen_3();
		if(licznik == 0) fprint(gen3,x3,N+3,3);
	// SFERA I KULA
		// sfera
 		double u1, u2, u3, u4, norm, skalar;
 		r = (double**)malloc(N*sizeof(double*));
 		for(int i=0; i<N; ++i) r[i] = (int*)malloc(3*sizeof(double));
 		x = (double*)malloc(N*sizeof(double));
 		y = (double*)malloc(N*sizeof(double));
 		z = (double*)malloc(N*sizeof(double));
		for(int i=0; i<N; ++i)
		{
			u1 = gen_3();
			u2 = gen_3();
			u3 = gen_3();
			u4 = gen_3();
			r[i][0] = x[i] = sqrt(-2*log(1-u1)) * cos(2*M_PI*u2);
			r[i][1] = y[i] = sqrt(-2*log(1-u1)) * sin(2*M_PI*u2);
			r[i][2] = z[i] = sqrt(-2*log(1-u3)) * cos(2*M_PI*u4);
				skalar = r[i][0]*r[i][0] + r[i][1]*r[i][1] + r[i][2]*r[i][2];
			norm = sqrt(skalar);
			r[i][0] = r[i][0]/norm;
			r[i][1] = r[i][1]/norm;
			r[i][2] = r[i][2]/norm;
			if(licznik == 0) fprintf(sfera, "%lf\t%lf\t%lf\n", r[i][0],r[i][1],r[i][2]);
		}
	
		// kula
		double ui, si;
		R = (double**)malloc(N*sizeof(double*));
 		for(int i=0; i<N; ++i) R[i] = (int*)malloc(3*sizeof(double));
		int d = 3;
		for(int i=0; i<N; ++i)
		{
			ui = gen_3();
			si = pow(ui,1.0/d);
			R[i][0] = r[i][0]*si;
			R[i][1] = r[i][1]*si;
			R[i][2] = r[i][2]*si;
			if(licznik == 0) fprintf(kula, "%lf\t%lf\t%lf\n", R[i][0],R[i][1],R[i][2]);
		}
	
	// GESTOSC KULI
		// dane
		int K = 10;
		double Delta = 1.0/K;
		int n[K+1];		// by indeksowac od 1
		int j;
		for(int i=0; i<=K; ++i) n[i]=0;
		// przedzialy
		for(int i=0; i<N; ++i)
		{
			skalar = R[i][0]*R[i][0] + R[i][1]*R[i][1] + R[i][2]*R[i][2];
			norm = sqrt(skalar);
			j = (int)(norm/Delta) + 1;
			++n[j];
		}
	/*	// sprawdzenie
		int sum=0;
		for(int i=1; i<=K; ++i) 
		{
			printf("%d\n", n[i]);
			sum+=n[i];
		}
	
		printf("sum = %d\n", sum);
	*/	// gestosc
		double rad[K+1], V[K+1], g[K+1];
		for(int i=0; i<=K; ++i) rad[i] = V[i] =  g[i] = 0;
	
		for(int i=1; i<=K; ++i)
		{
			rad[i] = Delta*i;
			V[i] = 4*M_PI/3*pow(rad[i],3);
			g[i] = n[i]/(V[i]-V[i-1]);
			fprintf(licznik==0 ? gest1 : licznik == 1 ? gest2 : gest3, "%d\t%d\t%lf\n", i, n[i], g[i]);
		}
		
		free(x1); free(x2); free(x3); free(x); free(y); free(z);
		for(int i=0;i<N;++i)
		{
			free(r[i]);
			free(R[i]);
		}
		free(r); free(R);
	}
		
	fclose(gen1);
	fclose(gen2);
	fclose(gen3);
	fclose(sfera);
	fclose(kula);
	fclose(gest1);
	fclose(gest2);
	fclose(gest3);
	return 0;
}

