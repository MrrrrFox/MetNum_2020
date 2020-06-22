#include <stdio.h>

#include "math.h"
#include <stdlib.h>
#include <time.h>

int main()
{
	srand(time(0));

	int N = 10000;
	float eps = 0.00001;
	int flaga = 0;
	int ile_iter = 0;
// stale zmienne:
	float beta = 0.4;
	float F_0 = 0.1;
	float Omega = 0.8;
// zmienne zmienne:
	float v_0 = 0.0;
	float h = 0.02;
	float omega = 1.0;

	float d_0[N];
	float d_1[N];
	float d_2[N];

	float x_n[N];
	float x_s[N];
	
	float b[N];

	float a_1 = 1.0;
	float a_2 = omega*omega*h*h - 2 - beta*h;
	float a_3 = 1 + beta*h;

	d_0[0] = d_0[1] = 1.0;
	d_1[0] = 0.0;
	d_1[1] = -1.0;
	d_2[0] = d_2[1] = 0.0;
	for(int i = 2; i<N; ++i)
	{
		d_0[i] = a_3;
		d_1[i] = a_2;
		d_2[i] = a_1;
	}

	b[0] = 1.0;
	for(int i = 1; i<N; ++i)
	{
		b[i] = F_0 * sin(Omega*h*i) * h * h;
	}
	//printf("d_0\t\td_1\t\td_2\t\tb\n");
	for(int i = 0; i<N; ++i)
	{
		//printf("%lf\t%lf\t%lf\t%lf\n", d_0[i], d_1[i], d_2[i], b[i]);
		x_s[i] = (float) rand() / RAND_MAX; 
		x_n[i] = 0.0;
	}

	while(ile_iter<100000)
	{
		for(int i = 0; i<N; ++i)
		{	
			if(i>=2)
			{
				x_n[i] = (1.0/d_0[i]) * (b[i] - d_1[i]*x_s[i-1] - d_2[i]*x_s[i-2]);
			}
			else
			{
				if(i == 1) x_n[i] = (1.0/d_0[i]) * (b[i] - d_1[i]*x_s[i-1]);
				else x_n[i] = (1.0/d_0[i]) * b[i];
			}
		}

		float sumas = 0.0;
		float suman = 0.0;
		
		for(int i = 0; i<N; ++i)
		{
			sumas += x_s[i] * x_s[i];
			suman += x_n[i] * x_n[i];
		}
		//printf("%lf\n", sumas-suman);
	
		if(fabs(suman - sumas) < eps) break;

		for(int i = 0; i<N; ++i)
		{
			//printf("%lf\t\t%lf\n", x_s[i], x_n[i]);
			x_s[i] = x_n[i];
		}
		++ile_iter;
	}

	printf("Ilosc iteracji: %d\n", ile_iter);

	//for(int i = 0; i<N; ++i)
	//{
	//	printf( "%lf\t%f\n", i*h, x_n[i] );
	//}

	return 0;
}

