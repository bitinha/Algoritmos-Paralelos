#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "travelingGreedy.c"

int main(int argc, char const *argv[])
{
	
	srand(0);

	// Number of towns to visit
	int N = 10;

	double *x = malloc(sizeof(double)*N);
	double *y = malloc(sizeof(double)*N);
	double *D = malloc(sizeof(double)*N*N);

	for (int i = 0; i < N; ++i)
	{
		x[i] = ((double)rand() / (double)RAND_MAX)*10;
		y[i] = ((double)rand() / (double)RAND_MAX)*10;
	}

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			D[i*N+j]=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
		}
	}


	int *order;

	double distance = greedy(D,N,&order);




	printf("%f\n", distance);

	printf("%d\n", order[1]);

	printf("%d\n", order[0]);

	free(order);
	free(D);
	free(y);
	free(x);

	return 0;
}