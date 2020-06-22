#include <stdio.h>
#include <complex.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////// WYPISYWANIE TABLICY
void print_complex(double complex * a,int rozmiar)
{
	for(int i = 0; i<rozmiar; ++i)
	{
		printf("%.1f\t%+.1fI\n", creal(a[i]), cimag(a[i]));
	}
}

////////////////////////////////////////////////////////////////////// LICZENIE Rj
double complex do_Rj(double complex * a, double complex * b, int rozmiar, double complex z0)
{
	b[rozmiar] = 0. + 0.I;
	for(int i = rozmiar-1; i>=0; --i)
	{
		b[i] = a[i+1] + z0 * b[i+1];
	}
	return a[0] + z0 * b[0];
}

////////////////////////////////////////////////////////////////////// DEFLACJA
void b_to_a(double complex * a,double complex * b, int rozmiar)
{
	for(int i = 0; i<rozmiar; ++i)
	{
		a[i] = b[i];
	}
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
	FILE *out1;
	out1 = fopen("out1.txt","w");
	FILE *out2;
	out2 = fopen("out2.txt","w");
	
	int n = 4;
	
	double complex *a;
	a = (double complex *) malloc( (n+1) * sizeof(double complex));
	a[0] = 16. + 8.I;
	a[1] = -20. + 14.I;
	a[2] = 4. - 8.I;
	a[3] = -4. + 1.I;
	a[4] = 1. + 0.I;
// TEST WYPISYWANIE
	//print_complex(a,n+1);

	double complex *b;
	b = (double complex *) malloc( (n+1) * sizeof(double complex));
	double complex *c;
	c = (double complex *) malloc( n * sizeof(double complex));
// TEST 1
	int stopien = n;
	double complex z0 = 0. + 0.I;
	double complex Rj = do_Rj(a,b,stopien,z0);
	//print_complex(b,n+1);
	//printf("Rj = %.1f %+.1fI\n\n", creal(Rj), cimag(Rj));
// TEST 2
	stopien = n-1;
	z0 = 0. + 0.I; // brak zmiany wartosci, ale by pokazac dla jakiego z0 to robimy
	double complex Rj_p = do_Rj(b,c,stopien,z0); // p jak prim
	//print_complex(c,n);
	//printf("Rj' = %.1f %+.1fI\n\n", creal(Rj_p), cimag(Rj_p));
	
	int IT_MAX = 20;
	stopien = n;
// TEST JEDEN PIERWIASTEK
	double complex z_first = z0;
	double complex z_next;
	for(int j = 1; j<= IT_MAX; ++j)
	{
		Rj = do_Rj(a,b,stopien,z_first);
		Rj_p = do_Rj(b,c,stopien-1,z_first);
		z_next = z_first - Rj/Rj_p;
		//printf("%d\t%.1f%+.1fI\t%.1f%+.1fI\t%.1f%+.1fI\n", j, creal(Rj), cimag(Rj), creal(Rj_p), cimag(Rj_p), creal(z_next), cimag(z_next));
		z_first = z_next;
	}
	//printf("zj = %.1f %+.1fI\n\n", creal(z_next), cimag(z_next));
	
// WLASCIWY ALGORYTM (wczesniejsze testy s¹ pozytywne)
// z0 = 0 + 0I
	z0 = 0. + 0.I;
	for(int i = n; i>=1; --i)
	{
		z_first = z0;
		fprintf(out1, "%.5f\t%.5f\n", creal(z_first), cimag(z_first));
		for(int j = 1; j<= IT_MAX; ++j)
		{
			Rj = do_Rj(a,b,i,z_first);
			Rj_p = do_Rj(b,c,i-1,z_first);
			z_next = z_first - Rj/Rj_p;
			fprintf(out1, "%.5f\t%.5f\n", creal(z_next), cimag(z_next));
			z_first = z_next;
		}
		fprintf(out1, "\n\n");
		b_to_a(a,b,n+1);
	}
// z0 = -10 - 10I
	a[0] = 16. + 8.I;
	a[1] = -20. + 14.I;
	a[2] = 4. - 8.I;
	a[3] = -4. + 1.I;
	a[4] = 1. + 0.I;
	z0 = -10. - 10.I;
	for(int i = n; i>=1; --i)
	{
		z_first = z0;
		fprintf(out2, "%.5f\t%.5f\n", creal(z_first), cimag(z_first));
		for(int j = 1; j<= IT_MAX; ++j)
		{
			Rj = do_Rj(a,b,i,z_first);
			Rj_p = do_Rj(b,c,i-1,z_first);
			z_next = z_first - Rj/Rj_p;
			fprintf(out2, "%.5f\t%.5f\n", creal(z_next), cimag(z_next));
			z_first = z_next;
		}
		fprintf(out2, "\n\n");
		b_to_a(a,b,n+1);
	}
	
	free(a);
	free(b);
	free(c);
	fclose(out1);
	fclose(out2);
	return 0;
}

