#include <stdio.h>

#include "math.h"

////////////////////////////////////////////////////////////////////// FUNKCJA 1
double fun1(double x)
{
	return log(pow(x,5)+3*pow(x,2)+x+9);
}

////////////////////////////////////////////////////////////////////// FUNKCJA 2
double fun2(double x)
{
	return pow(x,6);
}

////////////////////////////////////////////////////////////////////// ILOR ROZN 1
double F1(double (*fun)(double), double x1, double x2)
{
	return (fun(x2)-fun(x1))/(x2-x1);
}

////////////////////////////////////////////////////////////////////// ILOR ROZN 2
double F2(double (*fun)(double), double x1, double x2, double x3)
{
	return (F1(fun,x2,x3) - F1(fun,x1,x2))/(x3-x1);
}

////////////////////////////////////////////////////////////////////// MAIN
int main()
{
// dane
	double x_min = -1.5, x_max = 1., h = 0.01, minimum;
	double eps = 1e-7; // byl eps = 1e-10, ale powodowal komunikat #IND, #QNANO, wiec zmniejszylem dokladnosc do 1e-7 (przy 1e-8 bylo juz ok, 
						//ale zasugerowalem sie plikiem wyniki.pdf, gdzie mielismy wykres dla 7 iteracji dla fun1, a tyle dostajemy przy 10e-7)	 
						
////////////////////////////////////////////////////////////////////// wykres fun1 (log)
// plik
	FILE *logChart;
	logChart = fopen("logChart.txt","w");
// wykres
	for(double x = x_min; x<=x_max; x+=h)
	{
		fprintf(logChart, "%f\t%f\n", x, fun1(x));
	}
// zamkniecie pliku
	fclose(logChart);
////////////////////////////////////////////////////////////////////// fun1 dla x1 = -0.5
// punkty startowe - malejace wartosci funkcji
	printf("#################### FUN1, x1=-0.5 ####################\n");
	double x1 = -0.5; 
	double x2 = fun1(x1+h)<fun1(x1) ? x1+h : x1-h;
	double x3 = fun1(x2+h)<fun1(x2) ? x2+h : x2-h;
// plik
	FILE *fun1ad1;
	fun1ad1 = fopen("fun1ad1.txt","w");
// przybli¿enie x_m
	double x_m;
// dane c.d.
	int n = 10;
	
	for(int i = 0; i<n; ++i)
	{
		x_m = ((x1+x2)/2.0) - (F1(fun1,x1,x2)/(2*F2(fun1,x1,x2,x3)));
		fprintf(fun1ad1,"%i\t%f\t%f\t%f\t%f\t%f\t%f\n", i+1, x1,x2,x3,x_m,F1(fun1,x1,x2),F2(fun1,x1,x2,x3));
		// najblizszy x_n - sprawdzenie warunku
		if(fabs(x_m-x1) < eps || fabs(x_m-x2) < eps || fabs(x_m-x3) < eps)
		{
			minimum = x_m;
			break;
		}
		if(fabs(x_m-x1) < fabs(x_m-x2))
		{
			if(fabs(x_m-x2) < fabs(x_m-x3)) x3 = x_m;
			else x2 = x_m;
		}
		else
		{
			if(fabs(x_m-x1) < fabs(x_m-x3)) x3 = x_m;
			else x1 = x_m;
		}
	}
// zamkniecie pliku
	fclose(fun1ad1);
////////////////////////////////////////////////////////////////////// fun1 dla x1 = -0.9
	printf("\n\n#################### FUN1, x1=-0.9 ####################\n");
	x1 = -0.9; 
	x2 = fun1(x1+h)<fun1(x1) ? x1+h : x1-h;
	x3 = fun1(x2+h)<fun1(x2) ? x2+h : x2-h;
// plik
	FILE *fun1ad2;
	fun1ad2 = fopen("fun1ad2.txt","w");

	for(int i = 0; i<n; ++i)
	{
		x_m = ((x1+x2)/2.0) - (F1(fun1,x1,x2)/(2*F2(fun1,x1,x2,x3)));
		fprintf(fun1ad2,"%i\t%f\t%f\t%f\t%f\t%f\t%f\n", i+1, x1,x2,x3,x_m,F1(fun1,x1,x2),F2(fun1,x1,x2,x3));
		// najblizszy x_n - sprawdzenie warunku
		if(fabs(x_m-x1) < eps || fabs(x_m-x2) < eps || fabs(x_m-x3) < eps)
		{
			minimum = x_m;
			break;
		}
		if(fabs(x_m-x1) < fabs(x_m-x2))
		{
			if(fabs(x_m-x2) < fabs(x_m-x3)) x3 = x_m;
			else x2 = x_m;
		}
		else
		{
			if(fabs(x_m-x1) < fabs(x_m-x3)) x3 = x_m;
			else x1 = x_m;
		}
	}
// zamkniecie pliku
	fclose(fun1ad2);
////////////////////////////////////////////////////////////////////// fun2 dla x1 = 1.5
	printf("\n\n#################### FUN2, x1=1.5 ####################\n");
	x1 = 1.5; 
	x2 = fun1(x1+h)<fun1(x1) ? x1+h : x1-h;
	x3 = fun1(x2+h)<fun1(x2) ? x2+h : x2-h;
// plik
	FILE *fun2ad1;
	fun2ad1 = fopen("fun2ad1.txt","w");
// dane c.d.
	n = 100;

	for(int i = 0; i<n; ++i)
	{
		x_m = ((x1+x2)/2.0) - (F1(fun2,x1,x2)/(2*F2(fun2,x1,x2,x3)));
		fprintf(fun2ad1,"%i\t%f\t%f\t%f\t%f\t%f\t%f\n", i+1, x1,x2,x3,x_m,F1(fun2,x1,x2),F2(fun2,x1,x2,x3));
		// najblizszy x_n - sprawdzenie warunku
		if(fabs(x_m-x1) < eps || fabs(x_m-x2) < eps || fabs(x_m-x3) < eps)
		{
			minimum = x_m;
			break;
		}
		if(fabs(x_m-x1) < fabs(x_m-x2))
		{
			if(fabs(x_m-x2) < fabs(x_m-x3)) x3 = x_m;
			else x2 = x_m;
		}
		else
		{
			if(fabs(x_m-x1) < fabs(x_m-x3)) x3 = x_m;
			else x1 = x_m;
		}
	}
// zamkniecie pliku
	fclose(fun2ad1);

	return 0;
}

