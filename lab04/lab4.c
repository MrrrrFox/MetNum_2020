#include <stdio.h>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

#include "math.h"

double x_i(int i, double L, double dx)
{
	return -L/2. + dx * (i+1); 
}

double p_i(double alpha, double x_i)
{
	return 1+4*alpha*x_i*x_i;
}

void fill_Matrix(gsl_matrix* A, int n, double N, double dx, double alpha, double L)
{
	for(int i = 0; i<n; ++i)
	{
		gsl_matrix_set(A, i, i, 2*N/(dx*dx*p_i(alpha,x_i(i,L,dx))));
	}
	for(int i = 0; i<n-1; ++i)
	{
		gsl_matrix_set(A, i, i+1, -N/(dx*dx*p_i(alpha,x_i(i,L,dx))));
	}
	for(int i = 1; i<n; ++i)
	{
		gsl_matrix_set(A, i, i-1, -N/(dx*dx*p_i(alpha,x_i(i,L,dx))));
	}
}

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

int main()
{
	FILE *eigen_vals;
	FILE *alpha0;
	FILE *alpha100;
	
	eigen_vals = fopen("eigen_vals.txt","w");
	alpha0 = fopen("alpha0.txt","w");
	alpha100 = fopen("alpha100.txt","w");
	
	double L = 10;
	double N = 1;
	int n = 6;
	double dx = L/(n+1);
	double alpha;
//2.1
	//printf("dx = %.5lf\n", dx);
	for(int l1 = 0; l1<n; ++l1)
	{
		//printf("x_%d = %.5lf\n", l1, x_i(l1,L,dx));
	}
//2.2
	gsl_matrix* A = gsl_matrix_calloc(n,n);

	alpha = 0;	
	fill_Matrix(A, n, N, dx, alpha, L);
	//print_Matrix(A,n);

	alpha = 100;
	fill_Matrix(A, n, N, dx, alpha, L);
	//print_Matrix(A,n);
//2.3
	n = 200;
	dx = L/(n+1);
	alpha = 0;
	double dalpha = 2;
	gsl_vector_complex* eval = gsl_vector_complex_calloc(n);
	gsl_matrix_complex *evec = gsl_matrix_complex_calloc(n,n);
	gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
	
	for(alpha = 0;alpha<=100;alpha+=dalpha)
	{
		gsl_matrix_free(A);
		A = gsl_matrix_calloc(n,n);
		fill_Matrix(A,n,N,dx,alpha,L);

		gsl_eigen_nonsymmv(A,eval,evec,w);
		gsl_eigen_nonsymmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
		
		if(alpha == 0)
		{
			double tab[n][7];
			for(int i =0; i<n; ++i)
			{
				tab[i][0] = x_i(i,L,dx);
			}
			for(int numer = 0; numer<6; ++numer)
			{
				for(int i=0; i<n; ++i)
				{
					gsl_complex cv = gsl_matrix_complex_get(evec, i, numer);
					tab[i][numer+1] = GSL_REAL(cv);
				}
			}
			for(int i=0;i<n;++i)
			{
				for(int j=0; j<7; ++j)
				{
					fprintf(alpha0,"%.5lf\t",tab[i][j]);
				}
				fprintf(alpha0,"\n");
			}	
		}
		if(alpha == 100)
		{
			double tab[n][7];
			for(int i =0; i<n; ++i)
			{
				tab[i][0] = x_i(i,L,dx);
			}
			for(int numer = 0; numer<6; ++numer)
			{
				for(int i=0; i<n; ++i)
				{
					gsl_complex cv = gsl_matrix_complex_get(evec, i, numer);
					tab[i][numer+1] = GSL_REAL(cv);
				}
			}
			for(int i=0;i<n;++i)
			{
				for(int j=0; j<7; ++j)
				{
					fprintf(alpha100,"%.5lf\t",tab[i][j]);
				}
				fprintf(alpha100,"\n");
			}
		}
		
		fprintf(eigen_vals,"%.5lf\t", alpha);
		for(int i = 0; i<6; ++i)
		{
			gsl_complex cval = gsl_vector_complex_get(eval, i);
			double val = GSL_REAL(cval);
			fprintf(eigen_vals,"%.5lf\t", sqrt(val));
		}
		fprintf(eigen_vals,"\n");
	}

//2.4
	
//2.5

	fclose(eigen_vals);
	fclose(alpha0);
	fclose(alpha100);

	gsl_matrix_free(A); 
	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	gsl_eigen_nonsymmv_free(w);
	return 0;
}

