#include <stdio.h>
// zakladając, że w obecnym katalogu są pliki nr.h, nrutil.h, nrutil.c, gaussj.c
// #include "nr.h"
// #include "nrutil.h"

// alternatywnie, pracując na Taurusie można skorzystać z poniższych dyrektyw
#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"

#include <stdlib.h>
#include <string.h>

int main()
{
	int n = 7;
	float w = 1; // k/m = 1 => w^2 = 1 => w = 1
	float h = 0.1;
	float v_0 = 0;
	float a = 1; // w tresci zadania "A"

	float **vec = matrix(1, n, 1, 1);
	float **A = matrix(1, n, 1, n);

	// zerowanie macierzy:
	for (int i = 1; i < n+1; ++i)
	{
		for (int j = 1; j < n+1; ++j)
		{
			A[i][j] = 0;
		}
		vec[i][1] = 0;
	}

	// wypelnianie macierzy:
	for(int i = 1; i<n+1; ++i)
	{
		A[i][i] = 1;
	}

	A[2][1] = -1;
	float x = w*w*h*h - 2; // w^2 * h^2 - 2
	for(int i = 3; i<n+1; ++i)
	{
		A[i][i-1] = x;
	}

	for(int i = 3; i<n+1; ++i)
	{
		A[i][i-2] = 1;
	}

	// wypelnianie wektora:
	vec[1][1] = a;
	vec[2][1] = v_0 * h;
/*
	// wypisanie macierzy:
	 for (int i = 1; i < n+1; ++i)
	{
	 	for (int j = 1; j < n+1; ++j)
		{
	 		printf( "%.2lf\t", A[i][j] );
	 	}
	 	printf( "\n");
	 }
*/
	// rozwiazanie ukladu rownan
	gaussj(A, n, vec, 1);

	// wypisanie wyniku	- wywolanie programu w nastepujacy sposob: ./test > out.txt dzieki czemu dane od razu zostana wypisane do pliku
	for (int i = 1; i < n+1; ++i)
	{
		printf( "%lf\t%.6lf\n", (i-1)*h, vec[i][1] );
	}

	free_matrix(vec, 1, n, 1, 1);
	free_matrix(A, 1, n, 1, n);
	return 0;
}

