#include <stdio.h>
#include "math.h"

#define x_min 0.0
#define x_max M_PI

////////////////////////////////////////////////////////////////////// SILNIA
double silnia(int n)
{
	double result = 1;
	for(double i = 2; i<=n;++i) result*=i;
	//printf("silnia(%d) = %lf\n",n, result);
	return result;
}

////////////////////////////////////////////////////////////////////// I_j
long double I_j(int j, int k, int m)
{
	long double sum = 0;
	sum += pow(-1,j) * (pow(k*x_max, 2*j+m+2)) / (pow(k,m+1) * silnia(2*j+1) * (2*j+m+2));
	sum -= pow(-1,j) * (pow(k*x_min, 2*j+m+2)) / (pow(k,m+1) * silnia(2*j+1) * (2*j+m+2));
	//printf("j = %d:\tI_%d = %lf\n",j,j,sum);
	return sum;
}

////////////////////////////////////////////////////////////////////// f
double f(int index, int m, int k, double h)
{
	double x = x_min + h*index;
	//printf("f_%d = %lf\n", index, pow(x,m) * sin(k*x));
	return pow(x,m) * sin(k*x);
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
// obsluga plikow
	FILE *calka;
	calka = fopen("calka.txt","w");
	FILE *modul;
	modul = fopen("modul.txt","w");
// dane1
	double I;
	int m_tab[]={0,1,5}, k_tab[]={1,1,5};
	int m, k;
// dla roznych m i k
	for(int licznik1 = 0; licznik1<sizeof(m_tab)/sizeof(m_tab[0]); ++licznik1)
	{
		// m i k
		m = m_tab[licznik1];
		k = k_tab[licznik1];
		// obliczanie calki rozwijajac w Taylora
		for(int l=0; l<=30; ++l)
		{
			I = 0;
			for(int j=0; j<=l; ++j)
			{
				I += I_j(j,k,m);
			}
			fprintf(calka,"%d\t%lf\n",l,I);
		}
		//dane2
		int n_tab[] = {11,21,51,101,201};
		int n;
		double C;
		// dla roznych n
		for(int licznik2 = 0; licznik2<sizeof(n_tab)/sizeof(n_tab[0]); ++licznik2)
		{
			// dane2 cd
			n = n_tab[licznik2];
			int p = (n-1)/2;
			double h = (x_max-x_min)/(n-1);
			C=0;
			// obliczanie calki metoda Simpsona
			for(int j=1; j<=p; ++j)
			{
				C += f(2*j-2,m,k,h) + 4*f(2*j-1,m,k,h) + f(2*j,m,k,h);
			}
			C *= h/3;
			// porownanie modulu roznicy
			fprintf(modul, "%d\t%.10lf\n",n, C-I >= 0.0 ? C-I : -(C-I));
		}
		fprintf(calka, "\n\n");
		fprintf(modul, "\n\n");
	}
// obsluga plikow
	fclose(calka);
	fclose(modul);
	return 0;
}

