#include <stdio.h>

#include "sinft.c"
#include "four1.c"
#include "realft.c"
#include "nr.h"
#include "nrutil.c"
#include "nrutil.h"

#include "stdlib.h"
#include "time.h"


////////////////////////////////////////////////////////////////////// ZMIENNA LOSOWA
float X()
{
	return rand()/(RAND_MAX+1.0);
}

////////////////////////////////////////////////////////////////////// ZMIENNA LOSOWA IMITUJ¥CA SZUM
float a()
{
	float Y = X();
	int sign = (Y>0.5 ? 1 : -1);
	return 2*sign*X();
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
// losowosc
	srand(time(NULL));
// pliki
	FILE *file1;
	file1 = fopen("file1.txt","w");
	FILE *file2;
	file2 = fopen("file2.txt","w");
	FILE *file3;
	file3 = fopen("file3.txt","w");
	FILE *file4;
	file4 = fopen("file4.txt","w");
// zmienne
	int k_table[] = {10,8,6};
for(int iter = 0; iter<sizeof(k_table)/sizeof(k_table[0]); ++iter)
{
	int k = k_table[iter];
	int n = pow(2,k);
	float w = 2*2*M_PI/n;
// ZADANIE 1
	float y0[n];
	float y[n];
	for(int i = 0; i<n; ++i)
	{
		y0[i] = sin(w*(i+1)) + sin(2*w*(i+1)) + sin(3*w*(i+1));
		y[i] = y0[i] + a();
		if(iter == 0) fprintf(file1, "%d\t%f\t%f\n", (i+1), y[i], y0[i]);
	}
// ZADANIE 2
// transformata
	sinft(y,n);
	for(int i = 0; i<n; ++i)
	{
		if(iter == 0) fprintf(file2, "%d\t%f\n", (i+1), y[i]);
		if(i<=100 && iter==0) fprintf(file3, "%d\t%f\n", (i+1), y[i]);
	}
//ZADANIE 3
	float dysk = 0.25, max =-1e5;
// szukanie max
	for(int i = 0; i<n; ++i)
	{
		if(y[i]>max) max = y[i];
	}
// dyskryminacja
	dysk = dysk*max;
	for(int i = 0; i<n; ++i)
	{
		if(y[i]<dysk) y[i]=0;
	}
//ZADANIE 4
// transformata odwrotna
	sinft(y,n);
// przemnozenie
	for(int i = 0; i<n; ++i)
	{
		y[i] *= 2.0/n;
		fprintf(file4, "%d\t%f\t%f\n", (i+1), y[i], y0[i]);
	}
	fprintf(file4, "\n\n");
}
// zamkniecie plikow
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	return 0;
}

