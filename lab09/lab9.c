#include <stdio.h>

#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

#include "math.h"
#include "stdlib.h"
#include "time.h"

////////////////////////////////////////////////////////////////////// FUNKCJA g
float g(float x, float x_0, float sigma, float a0, float a1, float a2)
{
	// oryginalny wzor
	//return exp(-pow(x-x_0,2)/(2*pow(sigma,2)));
	
	// wielomianowy wzor
	return exp(a0 + a1*x+a2*pow(x,2));
}

////////////////////////////////////////////////////////////////////// FUNKCJA g2
float g2(float x, float x_0, float sigma, float a0, float a1, float a2)
{
	float X = rand()/(RAND_MAX+1.0);
	float alfa = 0.1;
	float delta = alfa *(X - 0.5);
	return g(x,x_0,sigma,a0,a1,a2)*(1+delta);
}

////////////////////////////////////////////////////////////////////// FUNKCJA G
float G_x(float x, float rozm, float ** b)
{
	float suma = 0;
	for(int i = 1; i<=rozm; ++i)
	{
		suma+= b[i][1]*pow(x,i-1);
	}
	return exp(suma);
}

////////////////////////////////////////////////////////////////////// FUNKCJA G2
float G2_x(float x, float rozm, float ** b)
{
	float suma = 0;
	for(int i = 1; i<=rozm; ++i)
	{
		suma+= b[i][1]*pow(x,i-1);
	}
	return exp(suma);
}

////////////////////////////////////////////////////////////////////// FUNKCJA f
float f(float x, float x_0, float sigma, float a0, float a1, float a2)
{
	return log(g(x,x_0,sigma,a0,a1,a2));
}

////////////////////////////////////////////////////////////////////// FUNKCJA f2
float f2(float x, float x_0, float sigma, float a0, float a1, float a2)
{
	return log(g2(x,x_0,sigma,a0,a1,a2));
}

////////////////////////////////////////////////////////////////////// WYPISYWANIE MACIERZY
void WypiszMacierz(float ** A, int rozm)
{
	printf("MACIERZ:\n");
	for(int i = 1; i<=rozm; ++i)
	{
		for(int j = 1; j<=rozm; ++j)
		{
			printf("%f\t", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

////////////////////////////////////////////////////////////////////// WYPISYWANIE WEKTORA
void WypiszWektor(float **vec, int rozm)
{
	printf("WEKTOR:\n");
	for(int i = 1; i<=rozm; ++i)
	{
		printf("%f\n", vec[i][1]);
	}
	printf("\n");
}

////////////////////////////////////////////////////////////////////// WYPISYWANIE WSPÓ£CZYNNIKÓW b_i
void Wypisz_b(float ** vec, int rozm)
{
	for(int i = 1; i<=rozm; ++i)
	{
		printf("b_%d = %.14f\t", i, vec[i][1]);
	}
	printf("\n");
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
	srand(time(NULL));
// zmienne uzywane w programie
	float x_0 = 2, sigma = 4, L = x_0 - 3*sigma, R = x_0 + 3* sigma;
	float a0 = -pow(x_0,2)/(2*pow(sigma,2)), a1 = x_0/pow(sigma,2), a2 = -1/(2*pow(sigma,2));
	int n_to_consider[] = {11,21,51,101};
	int n, m = 4;
// pliki
	FILE *first;
	FILE *second;
	first = fopen("fun_g.txt","w");
	second = fopen("apro_funG.txt","w");
	printf("#################### FUNKCJA g ####################\n");
	printf("a_0 = %f\t\ta_1 = %f\t\ta_2 = %f\n", a0,a1,a2);	
// DLA KAZDEGO n
	for(int l1 = 0; l1<sizeof(n_to_consider)/sizeof(n_to_consider[0]); ++l1)
	{
		n = n_to_consider[l1];
		if(n!=51) printf("n = %d\n", n);

		float h = (R-L)/(n-1);
// test1
		//printf("h = %f\n", h);
		for(int i = 1; i<=n; ++i)
		{
			float tmp = L+h*(i-1);
			//printf("x_%d = %f,\tf(x_%d) = %f\n", i, tmp, i, f(tmp,x_0,sigma,a0,a1,a2));
		}
		
// macierz i wektor
		float **r = matrix(1, m, 1, 1);
		float **G = matrix(1, m, 1, m);
// zerowanie macierzy i wektora
		for (int i = 1; i<=m; ++i)
		{
			for (int j = 1; j<=m; ++j)
			{
				G[i][j] = 0;
			}
			r[i][1] = 0;
		}
// wypelnianie macierzy
		float suma;
		for(int i = 1; i<=m; ++i)
		{
			for(int k = 1; k<=m; ++k)
			{
				suma = 0;
				for(float x = L; x<=R; x+=h)
				{
					suma += pow(x,i+k-2);
				}
				G[i][k] = suma;
			}
		}
// wypelnianie wektora
		for(int k = 1; k<=m; ++k)
		{
			suma = 0;
			for(float x = L; x<=R; x+=h)
			{
				suma += f(x,x_0,sigma,a0,a1,a2) * pow(x,k-1);
			}
			r[k][1] = suma;
		}
// test 2
		//WypiszMacierz(G,m);
		//WypiszWektor(r,m);	
// rozwiazanie ukladu rownan	
		gaussj(G, m, r, 1);
// wspolczynniki b
		if(n!=51) Wypisz_b(r,m);
// g vs G
		if(n==11 || n==101)
		{
			for(float x = L; x<=R; x+=0.01)
			{
				fprintf(first, "%f\t%f\n", x, g(x,x_0,sigma,a0,a1,a2));
			}
			fprintf(first, "\n\n");
			
			for(float x = L; x<=R; x+=0.01)
			{
				fprintf(second, "%f\t%f\n", x, G_x(x,m,r));
			}
			fprintf(second, "\n\n");	
		}
// obsluga pamieci
		free_matrix(r, 1, m, 1, 1);
		free_matrix(G, 1, m, 1, m);
	}
// zamkniecie plikow
	fclose(first);
	fclose(second);

/////////////////////////////////////////////////////////////////////////////////////////////////////// g2		
		
// pliki
	first = fopen("fun_g2.txt","w");
	second = fopen("apro_funG2.txt","w");	
	printf("\n\n#################### FUNKCJA g2 ####################\n");
	printf("a_0 = %f\t\ta_1 = %f\t\ta_2 = %f\n", a0,a1,a2);
// DLA KAZDEGO n
	for(int l1 = 0; l1<sizeof(n_to_consider)/sizeof(n_to_consider[0]); ++l1)
	{
		n = n_to_consider[l1];
		printf("n = %d\n", n);

		float h = (R-L)/(n-1);
// test1
		//printf("h = %f\n", h);
		for(int i = 1; i<=n; ++i)
		{
			float tmp = L+h*(i-1);
			//printf("x_%d = %f,\tf(x_%d) = %f\n", i, tmp, i, f(tmp,x_0,sigma,a0,a1,a2));
		}
		
// macierz i wektor
		float **r = matrix(1, m, 1, 1);
		float **G = matrix(1, m, 1, m);
// zerowanie macierzy i wektora
		for (int i = 1; i<=m; ++i)
		{
			for (int j = 1; j<=m; ++j)
			{
				G[i][j] = 0;
			}
			r[i][1] = 0;
		}
// wypelnianie macierzy
		float suma;
		for(int i = 1; i<=m; ++i)
		{
			for(int k = 1; k<=m; ++k)
			{
				suma = 0;
				for(float x = L; x<=R; x+=h)
				{
					suma += pow(x,i+k-2);
				}
				G[i][k] = suma;
			}
		}
// wypelnianie wektora
		for(int k = 1; k<=m; ++k)
		{
			suma = 0;
			for(float x = L; x<=R; x+=h)
			{
				suma += f2(x,x_0,sigma,a0,a1,a2) * pow(x,k-1);
			}
			r[k][1] = suma;
		}
// test 2
		//WypiszMacierz(G,m);
		//WypiszWektor(r,m);	
// rozwiazanie ukladu rownan	
		gaussj(G, m, r, 1);
// wspolczynniki b
		Wypisz_b(r,m);
// g vs G
		if(n==11 || n==101)
		{
			for(float x = L; x<=R; x+=h)
			{
				fprintf(first, "%f\t%f\n", x, g2(x,x_0,sigma,a0,a1,a2));
			}
			fprintf(first, "\n\n");
		
			for(float x = L; x<=R; x+=h)
			{
				fprintf(second, "%f\t%f\n", x, G2_x(x,m,r));
			}
			fprintf(second, "\n\n");
		}
// obsluga pamieci
		free_matrix(r, 1, m, 1, 1);
		free_matrix(G, 1, m, 1, m);
	}
// zamkniecie plikow
	fclose(first);
	fclose(second);

	return 0;
}

