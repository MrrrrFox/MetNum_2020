#include <stdio.h>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

#include "math.h"

////////////////////////////////////////////////////////////////////// WYSWIETLANIE MACIERZY
void print_Matrix(gsl_matrix* A, int n)
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

////////////////////////////////////////////////////////////////////// NOWE X_K_I1
void do_x_k_i1(gsl_vector * x_k_i1, gsl_matrix * W, gsl_vector * x_k_i0, int n)
{
	double suma;
	for(int i = 0; i<n; ++i)
	{
		suma = 0;
		for(int j = 0; j<n; ++j)
		{
			suma += gsl_matrix_get(W, i, j) * gsl_vector_get(x_k_i0, j);
		}
		gsl_vector_set(x_k_i1, i, suma);
	}
}

////////////////////////////////////////////////////////////////////// OBLICZANIE PRZYBLIZENIA LAMBDY
double do_lambda_k_i(gsl_vector * x_k_i1, gsl_vector * x_k_i0, int n)
{
	double licznik = 0;
	double mianownik = 0;
	for(int i = 0; i<n; ++i)
	{
		licznik += gsl_vector_get(x_k_i1, i) * gsl_vector_get(x_k_i0, i);
		mianownik += gsl_vector_get(x_k_i0, i) * gsl_vector_get(x_k_i0, i);
	}
	return licznik/mianownik;
}

////////////////////////////////////////////////////////////////////// NOWE X_K_I0
void do_x_k_i0(gsl_vector * x_k_i0, gsl_vector * x_k_i1, int n)
{
	double mianownik = 0;
	
	for(int i = 0; i<n; ++i)
	{
		mianownik += gsl_vector_get(x_k_i1, i) * gsl_vector_get(x_k_i1, i);
	}
	
	mianownik = sqrt(mianownik);

	for(int i = 0; i<n; ++i)
	{
		gsl_vector_set(x_k_i0, i, gsl_vector_get(x_k_i1, i) / mianownik);
	}
}

////////////////////////////////////////////////////////////////////// NOWA MACIERZ W
void update_W(gsl_matrix * W, double lambda_k_i, gsl_vector * x_k_i0, int n)
{
	gsl_matrix* W_copy = gsl_matrix_calloc(n,n);
	
	for(int i = 0; i<n; ++i)
	{
		for(int j = 0; j<n; ++j)
		{
			gsl_matrix_set(W_copy, i, j, gsl_matrix_get(W, i, j));
			
			gsl_matrix_set(W, i, j, gsl_matrix_get(W_copy, i, j) - lambda_k_i * gsl_vector_get(x_k_i0,i) * gsl_vector_get(x_k_i0,j) );
		}
	}
	
	gsl_matrix_free(W_copy); 
}

////////////////////////////////////////////////////////////////////// X_trans * A
void D1(gsl_matrix * D, gsl_matrix * X, gsl_matrix * A, int n)
{
	double suma;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<n;++j)
		{
			suma = 0;
			for(int k=0; k<n; ++k)
			{
				suma += gsl_matrix_get(X, k, i) * gsl_matrix_get(A, k, j);
			}
			gsl_matrix_set(D, i, j, suma);
		}
	}
}
	
////////////////////////////////////////////////////////////////////// D1 * X
void D2(gsl_matrix * D, gsl_matrix * X, int n)
{
	gsl_matrix* D_copy = gsl_matrix_calloc(n,n);
	
	for(int i = 0; i<n; ++i)
	{
		for(int j = 0; j<n; ++j)
		{
			gsl_matrix_set(D_copy, i, j, gsl_matrix_get(D, i, j));
		}
	}
	
	double suma;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<n;++j)
		{
			suma = 0;
			for(int k=0; k<n; ++k)
			{
				suma += gsl_matrix_get(D_copy, i, k) * gsl_matrix_get(X, k, j);
			}
			gsl_matrix_set(D, i, j, suma);
		}
	}
	
	gsl_matrix_free(D_copy); 
	
	
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
	FILE *matrixD;
	FILE *eigenval;
	
	matrixD = fopen("matrixD.txt","w");
	eigenval = fopen("eigenval.txt","w");
	

	int n = 7;

// 1.1 wypelnienie macierzy A
	gsl_matrix* A = gsl_matrix_calloc(n,n);
	for(int i=0;i<n;++i)
	{
		for(int j=0; j<n; ++j)
		{
			gsl_matrix_set(A, i, j, 1/(sqrt(2+abs(i-j))));
		}
	}
	//print_Matrix(A,n);
	
// 1.2 obliczanie wartosci wlasnych
	int IT_MAX = 12;
	int K_val = n;
	gsl_matrix* W = gsl_matrix_calloc(n,n);
	// W = A
	for(int i=0;i<n;++i)
	{
		for(int j=0; j<7; ++j)
		{
			gsl_matrix_set(W, i, j, 1/(sqrt(2+abs(i-j))));
		}
	}
	// inicjalizacja wektorów x, lambd oraz pomocniczej macierzy X
	gsl_vector * x_k_i0 = gsl_vector_calloc(n);
	gsl_vector * x_k_i1 = gsl_vector_calloc(n);
	double lambda_k_i;
	gsl_matrix* X = gsl_matrix_calloc(n,n);
	// algorytm
	for(int k = 0; k<K_val; ++k)
	{
		// poczatek iteracja dla kazdego k -> ustawianie poczatkowych wartosci wektora x
		for(int i = 0; i<n; ++i)
		{
			gsl_vector_set(x_k_i0, i, 1);
		}
		for(int i = 1; i<=IT_MAX; ++i)
		{
			// wyliczenie x_k_i1
			do_x_k_i1(x_k_i1, W, x_k_i0, n);
			// wyliczenie lambda_k_i
			lambda_k_i = do_lambda_k_i(x_k_i1, x_k_i0, n);
			// wyliczenie nowego x_k_i0
			do_x_k_i0(x_k_i0, x_k_i1, n);
			
			// wpisywanie lambd do pliku
			fprintf(eigenval, "%d\t%lf\n",i,lambda_k_i);
		}
		// aktualizacja W po obliczeniu jednej wartosci wlasnej
		update_W(W, lambda_k_i, x_k_i0, n);

		// formatowanie zapisu lambd w pliku
		fprintf(eigenval, "\n\n");
		
		// zapamietywanie wektorow wlasnych
		for(int i = 0; i<n; ++i)
		{
			gsl_matrix_set(X, i, k, gsl_vector_get(x_k_i0, i));
		}
	}

// 1.3 macierz D
	gsl_matrix* D = gsl_matrix_calloc(n,n);
	D1(D,X,A,n);
	D2(D,X,n);
	// wpisanie D do pliku
	for(int i=0;i<n;++i)
	{
		for(int j=0; j<7; ++j)
		{
			fprintf(matrixD,"%lf\t", gsl_matrix_get(D,i,j) );
		}
		fprintf(matrixD,"\n");
	}
	

	fclose(matrixD);
	fclose(eigenval);

	gsl_matrix_free(A); 
	gsl_matrix_free(W); 
	gsl_matrix_free(X); 
	gsl_matrix_free(D);
	gsl_vector_free(x_k_i0);
	gsl_vector_free(x_k_i1);

	return 0;
}

