#include <stdio.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"

#include "math.h"

#define x_min -5.
#define x_max 5.

////////////////////////////////////////////////////////////////////// FUNKCJA F1
float fun1(float x)
{
	return 1/(1+x*x);
}

////////////////////////////////////////////////////////////////////// FUNKCJA F2
float fun2(float x)
{
	return cos(2*x);
}

////////////////////////////////////////////////////////////////////// SKLEJKA KUBICZNA
float sklejka(int i, float x, float h)
{
	float x_i = x_min + (i-1)*h;
	if(x >= x_i-2*h && x < x_i-h) return 1/(h*h*h)*(x-(x_i-2*h))*(x-(x_i-2*h))*(x-(x_i-2*h));
	if(x >= x_i-h && x<x_i) return 1/(h*h*h)*(h*h*h + 3*h*h*(x-(x_i-h)) + 3*h*(x-(x_i-h))*(x-(x_i-h)) - 3*(x-(x_i-h))*(x-(x_i-h))*(x-(x_i-h)));
	if(x >= x_i && x<x_i+h) return 1/(h*h*h)*(h*h*h + 3*h*h*(x_i+h-x) + 3*h*(x_i+h-x)*(x_i+h-x) - 3*(x_i+h-x)*(x_i+h-x)*(x_i+h-x));
	if(x >= x_i+h && x<x_i+2*h) return 1/(h*h*h)*(x_i+2*h-x)*(x_i+2*h-x)*(x_i+2*h-x);
	if(x < x_min-3*h && x > x_max+3*h) return 0;
	return 0;
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
// zmienne uzywane w programie
	int n1[4] = {5,6,10,20};
	int n2[3] = {6,7,14};
	float eps = 0.01, h, tmp, suma;
	int n;
////////////////////////////////////////////////////////////////////// dla fun1
// pliki
	FILE *firstfun;
	FILE *firstpoints;
	firstfun = fopen("firstfun.txt","w");
	firstpoints = fopen("firstpoints.txt","w");
// dla kazdego wymaganego w zadaniu n
	for(int licznik1 = 0; licznik1<4; ++licznik1)
	{
		n = n1[licznik1];
// wyznaczenie h
		h = (x_max-x_min)/(n-1);
		//printf("h=%f\n",h);
// obliczenie x i f(x) oraz zapisanie wartosci do pliku
		for(int i = 0; i<n; ++i)
		{
			tmp = x_min+i*h;
			fprintf(firstpoints, "%f\t%f\n",tmp,fun1(tmp));
		}
		fprintf(firstpoints, "\n\n");
// macierz i wektory
		float **vec = matrix(1, n, 1, 1);
		float **vec2 = matrix(1, n+2, 1, 1);
		float **A = matrix(1, n, 1, n);
// zerowanie macierzy i wektora
		for (int i = 1; i < n+1; ++i)
		{
			for (int j = 1; j < n+1; ++j)
			{
				A[i][j] = 0;
			}
			vec[i][1] = 0;
		}
// wypelnianie macierzy
		for(int i = 1; i<n+1; ++i)
		{
			A[i][i] = 4;
		}
		A[1][2] = A[n][n-1] = 2;
		for(int i = 2; i<n; ++i)
		{
			A[i][i-1] = A[i][i+1] = 1;
		}
// pochodna
		float alfa = (fun1(x_min+eps) - fun1(x_min-eps))/(2*eps);
		float beta = (fun1(x_max+eps) - fun1(x_max-eps))/(2*eps);
		//printf("alfa=%f\tbeta=%f\n", alfa, beta);
// wypelnianie wektora
		vec[1][1] = fun1(x_min) + h/3*alfa;
		for(int i = 1; i<n-1; ++i)
		{
			vec[i+1][1] = fun1(x_min+i*h);
		}
		vec[n][1] = fun1(x_max) - h/3*beta;
	
		for (int i = 1; i < n+1; ++i)
		{
			//printf( "b_%d=%f\n", i, vec[i][1] );
		}
// rozwiazanie ukladu rownan	
		gaussj(A, n, vec, 1);
// wyznaczenie wartosci drugiego wektora
		for(int i = 1; i<n+1; ++i)
		{
			vec2[i+1][1]=vec[i][1];
		}
		vec2[1][1] = vec2[3][1] - h/3*alfa;
		vec2[n+2][1] = vec2[n][1] + h/3*beta;
		
		for (int i = 1; i < n+3; ++i)
		{
			//printf( "c_%d=%f\n", i-1, vec2[i][1] );
		}
// wreszcie interpolacja
		for(float x = x_min; x<=x_max; x+=eps)
		{
			suma = 0; 
// wyznaczanie sumy
			for(int i = 0; i<=n+1; ++i)
			{
				suma += vec2[i+1][1] * sklejka(i,x,h);
			}
			fprintf(firstfun, "%f\t%f\n", x, suma);
		}
		fprintf(firstfun, "\n\n");
// obsluga pamieci
		free_matrix(vec2, 1, n, 1, 1);
		free_matrix(vec, 1, n, 1, 1);
		free_matrix(A, 1, n, 1, n);
	}
// zamkniecie plikow
	fclose(firstfun);
	fclose(firstpoints);	
////////////////////////////////////////////////////////////////////// dla fun2
// obsluga plikow
	FILE *secondfun;
	FILE *secondpoints;
	secondfun = fopen("secondfun.txt","w");
	secondpoints = fopen("secondpoints.txt","w");
		
	for(int licznik1 = 0; licznik1<3; ++licznik1)
	{
		n = n2[licznik1];
// wyznaczenie h
		h = (x_max-x_min)/(n-1);
		//printf("h=%f\n",h);
// obliczenie x i f(x) oraz zapisanie wartosci do pliku
		for(int i = 0; i<n; ++i)
		{
			tmp = x_min+i*h;
			fprintf(secondpoints, "%f\t%f\n",tmp,fun2(tmp));
		}		
		fprintf(secondpoints, "\n\n");
// macierz i wektory
		float **vec = matrix(1, n, 1, 1);
		float **vec2 = matrix(1, n+2, 1, 1);
		float **A = matrix(1, n, 1, n);
// zerowanie macierzy i wektora
		for (int i = 1; i < n+1; ++i)
		{
			for (int j = 1; j < n+1; ++j)
			{
				A[i][j] = 0;
			}
			vec[i][1] = 0;
		}
// wypelnianie macierzy
		for(int i = 1; i<n+1; ++i)
		{
			A[i][i] = 4;
		}
		A[1][2] = A[n][n-1] = 2;
		for(int i = 2; i<n; ++i)
		{
			A[i][i-1] = A[i][i+1] = 1;
		}
// pochodna
		float alfa = (fun2(x_min+eps) - fun2(x_min-eps))/(2*eps);
		float beta = (fun2(x_max+eps) - fun2(x_max-eps))/(2*eps);
		//printf("alfa=%f\tbeta=%f\n", alfa, beta);
		
// wypelnianie wektora
		vec[1][1] = fun2(x_min) + h/3*alfa;
		for(int i = 1; i<n-1; ++i)
		{
			vec[i+1][1] = fun2(x_min+i*h);
		}
		vec[n][1] = fun2(x_max) - h/3*beta;
	
		for (int i = 1; i < n+1; ++i)
		{
			//printf( "b_%d=%f\n", i, vec[i][1] );
		}
// rozwiazanie ukladu rownan
		gaussj(A, n, vec, 1);
// wyznaczenie wartosci drugiego wektora
		for(int i = 1; i<n+1; ++i)
		{
			vec2[i+1][1]=vec[i][1];
		}
		vec2[1][1] = vec2[3][1] - h/3*alfa;
		vec2[n+2][1] = vec2[n][1] + h/3*beta;

		for (int i = 1; i < n+3; ++i)
		{
			//printf( "c_%d=%f\n", i-1, vec2[i][1] );
		}
// wreszcie interpolacja
		for(float x = x_min; x<=x_max; x+=eps)
		{
			suma = 0; 
// wyznaczanie sumy
			for(int i = 0; i<=n+1; ++i)
			{
				suma += vec2[i+1][1] * sklejka(i,x,h);
			}
			fprintf(secondfun, "%f\t%f\n", x, suma);
		}
		fprintf(secondfun, "\n\n");
// obsluga pamieci
		free_matrix(vec2, 1, n, 1, 1);
		free_matrix(vec, 1, n, 1, 1);
		free_matrix(A, 1, n, 1, n);
	}
// zamkniecie plikow
	fclose(secondfun);
	fclose(secondpoints);
	
	return 0;
}

