#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <mpi.h>

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

void plotGraph(char* alg, int N, int *order, double *x, double *y) {
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    int i;
    
    fprintf(gnuplotPipe, "set title 'travelingSalesman'\n");
    fprintf(gnuplotPipe, "set terminal svg\n");
    fprintf(gnuplotPipe, "set output '%s.svg'\n", alg);
    fprintf(gnuplotPipe, "plot '-' title 'Greedy' with linespoints\n");
    for (i = 0; i < N - 1; i++)
        fprintf(gnuplotPipe, "%g %g\n", x[order[i]], y[order[i]]);
    fprintf(gnuplotPipe, "%g %g\n", x[order[0]], y[order[0]]);
    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);
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


    int rank, n_proc;

    MPI_Status status;
    MPI_Init(&argc, &argv); //Initializes the library
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get the id of the current process


	srand(seed+rank);

	int *order,*orderMC,*orderSA;

	int starting_town = rand() % N;

	double distance = greedy(D,N,&order,starting_town);




    
    
    double distanceMC = MC(D, N, &orderMC, iterations);


    double distanceSA = SA(D, N, &orderSA, iterations, 1);
	

    if(rank==0){
    	printf("--------------------------------------------\n");
    	printf("\n\nRESULTADO DO PROCESSO %d\n", rank);
		
		if (printPath)
		{
			printf("Distance Greedy: %f\n", distance);
        	printPathOutput(N, order, D);
            char *greedy = "Greedy";
            plotGraph(greedy, N, order, x, y);
			    
		    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
        	printPathOutput(N, orderMC, D);
             char *mc = "Monte-Carlo";
            plotGraph(mc, N, orderMC, x, y);
			
		    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
        	printPathOutput(N, orderSA, D);
            char *sa = "Simulated-Annealing";
            plotGraph(sa, N, orderSA, x, y);
			printf("\n\n");
			

	    	for (int i = 1; i < n_proc; ++i)
	    	{
	    		MPI_Recv(order,N,MPI_INT,i,0,MPI_COMM_WORLD,&status);
	    		MPI_Recv(orderMC,N,MPI_INT,i,1,MPI_COMM_WORLD,&status);
	    		MPI_Recv(orderSA,N,MPI_INT,i,2,MPI_COMM_WORLD,&status);

	    
	    	printf("--------------------------------------------\n");
    			printf("\n\nRESULTADO DO PROCESSO %d\n", i);

				printf("Distance Greedy: %f\n", distance);
	        	printPathOutput(N, order, D);
				    
			    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
	        	printPathOutput(N, orderMC, D);
				
			    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
	        	printPathOutput(N, orderSA, D);
			    printf("\n\n");
	    	}
	    }else{

			printf("Distance Greedy: %f\n", distance);
			    
		    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
			
		    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
			printf("\n\n");
			

	    	for (int i = 1; i < n_proc; ++i)
	    	{
	    		MPI_Recv(&distance,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&distanceMC,1,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&distanceSA,1,MPI_DOUBLE,i,2,MPI_COMM_WORLD,&status);
	    		MPI_Recv(order,N,MPI_INT,i,3,MPI_COMM_WORLD,&status);
	    		MPI_Recv(orderMC,N,MPI_INT,i,4,MPI_COMM_WORLD,&status);
	    		MPI_Recv(orderSA,N,MPI_INT,i,5,MPI_COMM_WORLD,&status);

	    
    			printf("--------------------------------------------\n");
    			printf("\n\nRESULTADO DO PROCESSO %d\n", i);

				printf("Distance Greedy: %f\n", distance);
				    
			    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
				
			    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
			    printf("\n\n");
	    	}

	    }
    }else{
    	MPI_Send(&distance,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		MPI_Send(&distanceMC,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
		MPI_Send(&distanceSA,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    	MPI_Send(order,N,MPI_INT,0,3,MPI_COMM_WORLD);
		MPI_Send(orderMC,N,MPI_INT,0,4,MPI_COMM_WORLD);
		MPI_Send(orderSA,N,MPI_INT,0,5,MPI_COMM_WORLD);
    }



	free(order);
	free(orderSA);
	free(orderMC);
	free(D);
	free(y);
	free(x);

    MPI_Finalize();


	return 0;
}
