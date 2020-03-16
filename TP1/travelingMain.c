#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "travelingGreedy.c"
#include "travelingMC.c"

int main(int argc, char const *argv[])
{
	int printPath = 0;
	int iterations = 100;
	int seed = 0;

	// Number of towns to visit
	int N = 10;

	if (argc>1)
	{
		N = atoi(argv[1]);
	}

	if (argc>2)
	{
		iterations = atoi(argv[2]);
	}

	if (argc>3)
	{
		seed = atoi(argv[3]);
	}
	
	srand(seed);

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

	int starting_town = rand() % N;

	double distance = greedy(D,N,&order,starting_town);




	printf("distance Greedy: %f\n", distance);

	if(printPath){
		for (int i = 0; i < N-1; ++i)
		{
			printf("Town %d to town %d: %f\n", order[i], order[i+1], D[order[i]*N+order[i+1]]);
			
		}
		printf("Town %d to town %d: %f\n", order[N-1], order[0], D[order[N-1]*N+order[0]]);
	}
    
    
    double distanceMC = MC(D, N, &order, iterations);

    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
    
    if(printPath){
        for (int i = 0; i < N-1; ++i) {
    		printf("Town %d to town %d: %f\n", order[i], order[i+1], D[order[i]*N+order[i+1]]);
    	}
    	printf("Town %d to town %d: %f\n", order[N-1], order[0], D[order[N-1]*N+order[0]]);
    }


    double distanceSA = SA(D, N, &order, iterations, 1);
    
    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
    if(printPath){
        for (int i = 0; i < N-1; ++i) {
    		printf("Town %d to town %d: %f\n", order[i], order[i+1], D[order[i]*N+order[i+1]]);
    	}
    	printf("Town %d to town %d: %f\n", order[N-1], order[0], D[order[N-1]*N+order[0]]);
    }
	
	free(order);
	free(D);
	free(y);
	free(x);

	return 0;
}
