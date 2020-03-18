#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>

#include "travelingGreedy.c"
#include "travelingMC.c"

void usage(const char *path) {
    // take only the last portion of the path
    const char *basename = strrchr ( path, '/' );
    basename = basename ? basename + 1 : path;

    printf("usage: %s [OPTION]\n", basename);
    printf(" -h, --help\t\t\tPrint this help and exit.\n");
    printf(" -n, --number-of-towns NUMBER\tSelect the number of towns to visit.\n");
    printf(" -i, --iterations NUMBER\tSelect the number of iterations.\n");
    printf(" -s, --seed SEED\t\tUse a custom seed.\n");
    printf(" -p, --printPath]\t\tPrint path of all algorithms.\n");
}

void printPathOutput(int N, int *order, double *D) {
    for (int i = 0; i < N-1; ++i) {
        printf("Town %d to town %d: %f\n", order[i], order[i+1], D[order[i]*N+order[i+1]]);
    }
    printf("Town %d to town %d: %f\n", order[N-1], order[0], D[order[N-1]*N+order[0]]);
}

int main(int argc, char *argv[])
{
	int printPath = 0;
	int iterations = 100;
	int seed = 0;

	// Number of towns to visit
	int N = 10;

	struct option longopts[] = {
        { "help", no_argument, NULL, 'h' },
        { "number-of-towns", required_argument, NULL, 'n' },
        { "iterations", required_argument, NULL, 'i' },
        { "seed", required_argument, NULL, 's' },
        { "printPath", no_argument, NULL, 'p' }
    };

    while (1) {
        int opt = getopt_long ( argc, argv, "hpn:i:s:", longopts, 0 );
        if (opt == -1) {
            /* a return value of -1 indicates that there are no more options */
            break;
        }
        switch (opt) {
            case 'h':
                usage(argv[0]);
                return 0;
            case 'n':
                N = atoi(optarg);
                break;
            case 'i':
                iterations = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'p':
                printPath = 1;
                break;
            case '?':
                usage(argv[0]);
                return 1;
            default:
                break;
        }
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
		printPathOutput(N, order, D);
	}
    
    
    double distanceMC = MC(D, N, &order, iterations);

    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
    
    if(printPath){
        printPathOutput(N, order, D);
    }


    double distanceSA = SA(D, N, &order, iterations, 1);
    
    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
    if(printPath){
        printPathOutput(N, order, D);
    }
	
	free(order);
	free(D);
	free(y);
	free(x);

	return 0;
}
