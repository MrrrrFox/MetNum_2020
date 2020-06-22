#include <stdio.h>

#include "gsl/gsl_matrix.h"

#include "math.h"

////////////////////////////////////////////////////////////////////// WARTOSC FUNKCJI
double func(double x)
{
	return 1/(1+x*x);
}

////////////////////////////////////////////////////////////////////// WYSWIETLANIE MACIERZY
void print_Matrix(gsl_matrix * A, int n)
{
	for(int i = 0; i<n; ++i)
	{
		for(int j = 0; j<n; ++j)
		{
			printf( "%e\t", gsl_matrix_get(A, i, j) );
		}
		printf("\n");
	}
	printf("\n");
}

////////////////////////////////////////////////////////////////////// X_I DLA WIELOMIANU CZEBYSZEWA
double czebyszew(int n, double a, double b, int i)
{
	return 0.5 * ((a-b) * cos(M_PI*(2*i+1)/(2*n+2)) + (a+b));
}

////////////////////////////////////////////////////////////////////// ILORAZ ROZNICOWY
double ilor(gsl_matrix * A, int i, int j, double a, double dx)
{
	return (gsl_matrix_get(A,i,j-1)-gsl_matrix_get(A,i-1,j-1))/((a+i*dx)-(a+(i-j)*dx));
}

////////////////////////////////////////////////////////////////////// ILORAZ ROZNICOWY DLA WIELOMIANU CZEBYSZEWA
double ilor2(gsl_matrix * A, int n, int i, int j, double a, double b)
{
	return (gsl_matrix_get(A,i,j-1)-gsl_matrix_get(A,i-1,j-1)) / (czebyszew(n,a,b,i) - czebyszew(n,a,b,i-j));
}

////////////////////////////////////////////////////////////////////// FUNCKJA POMOCNICZA - ILOCZYN (X-X_I) DLA W_N
double iloczyn(double x, int j, double a, double dx)
{
	double iloczyn = 1;
	for(int i = 0; i<=j-1; ++i)
	{
		iloczyn *= (x-(a+i*dx)); 
	}
	return iloczyn;
}

////////////////////////////////////////////////////////////////////// FUNCKJA POMOCNICZA - ILOCZYN (X-X_I) DLA W_N2
double iloczyn2(double n, double x, int j, double a, double b)
{
	double iloczyn = 1;
	for(int i = 0; i<=j-1; ++i)
	{
		iloczyn *= (x-czebyszew(n,a,b,i)); 
	}
	return iloczyn;
}

////////////////////////////////////////////////////////////////////// WIELOMIAN INTERPOLUJACY
double Wn(gsl_matrix * A, int n,  double x, double a, double dx)
{
	double suma = 0;
	for(int j = 0; j<=n; ++j)
	{
		suma += gsl_matrix_get(A,j,j) * iloczyn(x,j,a,dx);
	}
	return suma;
}

////////////////////////////////////////////////////////////////////// WIELOMIAN INTERPOLUJACY (WERSJA DLA CZEBYSZEWA)
double Wn2(gsl_matrix * A, int n,  double x, double a, double b)
{
	double suma = 0;
	for(int j = 0; j<=n; ++j)
	{
		suma += gsl_matrix_get(A,j,j) * iloczyn2(n,x,j,a,b);
	}
	return suma;
}



////////////////////////////////////////////////////////////////////// MAIN
int main()
{
	double a = -5, b=5;
// wezly rownoodlegle
	FILE *normal;
	FILE *wezly1;
	normal = fopen("normal.txt","w");
	wezly1 = fopen("p1.txt","w");
	for(int n = 5; n<=20; n+=5)
	{
		// obliczenie dx
		double dx = (b-a)/n;
		gsl_matrix* A = gsl_matrix_calloc(n+1,n+1);
		
		// ustawienie kolumny o indeksie 0 i zapisanie wezlow
		for(int i = 0; i<=n; ++i)
		{
			gsl_matrix_set(A, i, 0, func(a+i*dx));
			fprintf(wezly1, "%lf\t%lf\n", a+i*dx, func(a+i*dx));
		}

		fprintf(wezly1, "\n\n");
	
		// ilorazy roznicowe
		for(int j = 1; j<=n;++j)
		{
			for(int i = j; i<=n; ++i)
			{
				gsl_matrix_set(A,i,j, ilor(A,i,j,a,dx));
			}
		}
		
		//print_Matrix(A,n+1);
		
		//wielomina interpolujacy
		for(double x = a; x<=b; x+=0.02)
		{
			fprintf(normal, "%lf\t%lf\n",x,Wn(A,n,x,a,dx));
		}
		
		fprintf(normal, "\n\n");
		
		gsl_matrix_free(A);
	}
	fclose(normal);
	fclose(wezly1);
	
// wezly zoptymalizowane
	FILE *opt;
	FILE *wezly2;
	opt = fopen("opt.txt","w");
	wezly2 = fopen("p2.txt", "w");
	for(int n = 5; n<=20; n+=5)
	{
		gsl_matrix* A = gsl_matrix_calloc(n+1,n+1);
		
		// ustawienie kolumny o indeksie 0 i zapisanie wezlow
		for(int i = 0; i<=n; ++i)
		{
			gsl_matrix_set(A, i, 0, func(czebyszew(n,a,b,i)));
			fprintf(wezly2, "%lf\t%lf\n",czebyszew(n,a,b,i), func(czebyszew(n,a,b,i)));
		}

		fprintf(wezly2, "\n\n");

		// ilorazy roznicowe
		for(int j = 1; j<=n;++j)
		{
			for(int i = j; i<=n; ++i)
			{
				gsl_matrix_set(A,i,j, ilor2(A,n,i,j,a,b));
			}
		}
		
		//print_Matrix(A,n+1);
		
		// wielomian interpolujacy
		for(double x = a; x<=b; x+=0.02)
		{
			fprintf(opt, "%lf\t%lf\n",x,Wn2(A,n,x,a,b));
		}
		
		fprintf(opt, "\n\n");
		
		gsl_matrix_free(A);
	}
	fclose(opt);
	fclose(wezly2);

	return 0;
}

