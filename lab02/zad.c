#include <stdio.h>

#include "gsl/gsl_linalg.h"

#include <stdlib.h>
#include <string.h>

int main()
{
	int N = 4;	// rozmiar macierzy
	float delta = 2;

	// deklaracja macierzy A, wektora permutacji p oraz zmiennej znasku sigma
	gsl_matrix * A = gsl_matrix_calloc(N, N);
	gsl_matrix * A2 = gsl_matrix_calloc(N, N);
	gsl_permutation * p = gsl_permutation_calloc(N);
	int signum;
	
	// macierz A
	for(int i = 0; i<N; ++i)
	{
		for(int j = 0; j<N; ++j)
		{
			gsl_matrix_set(A, i, j, 1.0/(i+j+delta));
			gsl_matrix_set(A2, i, j, 1.0/(i+j+delta));
			//printf( "%lf\t", gsl_matrix_get(A, i, j) );
		}
			//printf("\n");
	}

	double normaA = gsl_matrix_get(A, 0, 0);
	for(int i = 0; i<N; ++i)
	{
		for(int j = 0; j<N; ++j)
		{
			if(normaA < gsl_matrix_get(A, i, j)) 
			{
				normaA = gsl_matrix_get(A, i, j);
			}
		}
	}

	// dekompozycja
	gsl_linalg_LU_decomp(A, p, &signum);

	// deklaracja wektora b i x
	gsl_vector * b = gsl_vector_calloc(N);
	gsl_vector * x = gsl_vector_calloc(N);;

	// wyliczenie wyznacznika A
	double detA = 1;
	printf("Diagonala: ");
	for(int i = 0; i<N; ++i)
	{
		detA *= gsl_matrix_get(A, i, i);
		printf("%lf\t", gsl_matrix_get(A, i, i));
	}
	printf("\n\nWyznacznik A to %e\n\n", detA*signum);

	gsl_matrix * A_odw = gsl_matrix_calloc(N, N);
	// wyliczenie A^-1
	printf("A^-1:\n");
	for(int i = 0; i<N; ++i)
	{
		for(int j = 0; j<N; ++j)
		{
			gsl_vector_set(b, j, 0);
		}
		gsl_vector_set(b, i, 1);
		gsl_linalg_LU_solve(A, p, b, x);
		for(int j = 0; j<N; ++j)
		{
			gsl_matrix_set(A_odw, i, j, gsl_vector_get(x, j));
			printf( "%lf\t", gsl_matrix_get(A_odw, i, j) );
		}
		printf("\n");
	}
	
	printf("\n");


	gsl_matrix * iloczyn = gsl_matrix_calloc(N, N);
	// iloczyn macierzy - blad
	printf("A * A^-1:\n");
	double suma;
	for(int i = 0; i<N; ++i)
	{
		for(int j = 0; j<N; ++j)
		{
			suma = 0;
			for(int k = 0; k<N; ++k)
			{
				suma += gsl_matrix_get(A2, i, k) * gsl_matrix_get(A_odw, k, j);
			}
			gsl_matrix_set(iloczyn, i, j, suma);
		}
	}

	// wypisanie iloczynu
	for(int i = 0; i<N; ++i)
	{
		for(int j = 0; j<N; ++j)
		{
			printf( "%e\t", gsl_matrix_get(iloczyn, i, j) );
		}
		printf("\n");
	}
		
	// wskaznik uwarunkowania
	double normaA_odw = gsl_matrix_get(A_odw, 0, 0);
	for(int i = 0; i<N; ++i)
	{
		for(int j = 0; j<N; ++j)
		{
			if(normaA_odw < gsl_matrix_get(A_odw, i, j)) 
			{
				normaA_odw = gsl_matrix_get(A_odw, i, j);
			}
		}
	}
	double wsk = normaA * normaA_odw;
	printf("\nWskaznik uwarunkowania = %lf\n", wsk);

	// zwalnianie pamieci
	gsl_matrix_free(A); 
	gsl_matrix_free(A2); 
	gsl_matrix_free(A_odw);
	gsl_matrix_free(iloczyn);
	gsl_permutation_free(p); 
	gsl_vector_free(b);
	gsl_vector_free(x);  
	return 0;
}

