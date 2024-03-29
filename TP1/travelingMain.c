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
    fprintf(gnuplotPipe, "set yrange[0:10]\n");
    fprintf(gnuplotPipe, "set xrange[0:10]\n");
    fprintf(gnuplotPipe, "set output '%s.svg'\n", alg);
    fprintf(gnuplotPipe, "plot '-' title '%s' with linespoints, '-' title 'Starting Town' pt 5\n",alg);
    for (i = 0; i < N ; i++)
        fprintf(gnuplotPipe, "%g %g\n", x[order[i]], y[order[i]]);
    fprintf(gnuplotPipe, "%g %g\n", x[order[0]], y[order[0]]);
    fprintf(gnuplotPipe, "e\n");
    fprintf(gnuplotPipe, "%g %g\n",x[order[0]], y[order[0]]);
    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);
}

void printPathOutput(int N, int *order, double *x, double *y) {
    for (int i = 0; i < N-1; ++i) {
        printf("Town %d to town %d: %f\n", order[i], order[i+1], sqrt(pow(x[order[i]]-x[order[i+1]],2)+pow(y[order[i]]-y[order[i+1]],2)));
    }
    printf("Town %d to town %d: %f\n", order[N-1], order[0], sqrt(pow(x[order[N-1]]-x[order[0]],2)+pow(y[order[N-1]]-y[order[0]],2)));
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
/*
	double t_greedy = MPI_Wtime();
	double distance = greedyP(x,y,N,&order,starting_town);
    t_greedy = MPI_Wtime() - t_greedy;
    
	double t_MC = MPI_Wtime();
    double distanceMC = MCP(x, y, N, &orderMC, iterations);
    t_MC = MPI_Wtime() - t_MC;

	double t_SA = MPI_Wtime();
    double distanceSA = SAP(x, y, N, &orderSA, iterations, 100000);
	t_SA = MPI_Wtime() - t_SA;
	*/
	double t_greedy = MPI_Wtime();
	double distance = greedy(D,N,&order,starting_town);
    t_greedy = MPI_Wtime() - t_greedy;
    
	double t_MC = MPI_Wtime();
    double distanceMC = MC(D, N, &orderMC, iterations);
    t_MC = MPI_Wtime() - t_MC;

	double t_SA = MPI_Wtime();
    double distanceSA = SA(D, N, &orderSA, iterations, 1);
	t_SA = MPI_Wtime() - t_SA;

    if(rank==0){
    	printf("--------------------------------------------\n");
    	printf("\n\nRESULTADO DO PROCESSO %d\n", rank);
		
		if (printPath)
		{
			printf("Distance Greedy: %f\n", distance);
        	printPathOutput(N, order, x, y);
            char *greedy = "Greedy";
            plotGraph(greedy, N, order, x, y);
            printf("Tempo: %f\n", t_greedy);
			    
		    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
        	printPathOutput(N, orderMC, x, y);
            char *mc = "Monte-Carlo";
            plotGraph(mc, N, orderMC, x, y);
            printf("Tempo: %f\n", t_MC);

		    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
        	printPathOutput(N, orderSA, x, y);
            char *sa = "Simulated-Annealing";
            plotGraph(sa, N, orderSA, x, y);
            printf("Tempo: %f\n", t_SA);
			printf("\n\n");
			

	    	for (int i = 1; i < n_proc; ++i)
	    	{
	    		MPI_Recv(&distance,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&distanceMC,1,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&distanceSA,1,MPI_DOUBLE,i,2,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&t_greedy,1,MPI_DOUBLE,i,3,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&t_MC,1,MPI_DOUBLE,i,4,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&t_SA,1,MPI_DOUBLE,i,5,MPI_COMM_WORLD,&status);
	    		MPI_Recv(order,N,MPI_INT,i,6,MPI_COMM_WORLD,&status);
	    		MPI_Recv(orderMC,N,MPI_INT,i,7,MPI_COMM_WORLD,&status);
	    		MPI_Recv(orderSA,N,MPI_INT,i,8,MPI_COMM_WORLD,&status);

	    
	    	printf("--------------------------------------------\n");
    			printf("\n\nRESULTADO DO PROCESSO %d\n", i);

				printf("Distance Greedy: %f\n", distance);
	        	printPathOutput(N, order, x, y);
            	printf("Tempo: %f\n", t_greedy);
				    
			    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
	        	printPathOutput(N, orderMC, x, y);
            	printf("Tempo: %f\n", t_MC);
				
			    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
	        	printPathOutput(N, orderSA, x, y);
            	printf("Tempo: %f\n", t_SA);
			    printf("\n\n");
	    	}
	    }else{

			printf("Distance Greedy: %f\n", distance);
            printf("Tempo: %f\n", t_greedy);
			    
		    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
            printf("Tempo: %f\n", t_MC);
			
		    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
            printf("Tempo: %f\n", t_SA);
			printf("\n\n");
			

	    	for (int i = 1; i < n_proc; ++i)
	    	{
	    		MPI_Recv(&distance,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&distanceMC,1,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&distanceSA,1,MPI_DOUBLE,i,2,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&t_greedy,1,MPI_DOUBLE,i,3,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&t_MC,1,MPI_DOUBLE,i,4,MPI_COMM_WORLD,&status);
	    		MPI_Recv(&t_SA,1,MPI_DOUBLE,i,5,MPI_COMM_WORLD,&status);

	    
    			printf("--------------------------------------------\n");
    			printf("\n\nRESULTADO DO PROCESSO %d\n", i);

				printf("Distance Greedy: %f\n", distance);
            	printf("Tempo: %f\n", t_greedy);
				    
			    printf("\nDistance Monte-Carlo: %f\n", distanceMC);
            	printf("Tempo: %f\n", t_MC);
				
			    printf("\nDistance Simulated Annealing: %f\n", distanceSA);
            	printf("Tempo: %f\n", t_SA);
			    printf("\n\n");
	    	}

	    }
    }else{
    	MPI_Send(&distance,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		MPI_Send(&distanceMC,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
		MPI_Send(&distanceSA,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
    	MPI_Send(&t_greedy,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
		MPI_Send(&t_MC,1,MPI_DOUBLE,0,4,MPI_COMM_WORLD);
		MPI_Send(&t_SA,1,MPI_DOUBLE,0,5,MPI_COMM_WORLD);
		if(printPath){
	    	MPI_Send(order,N,MPI_INT,0,6,MPI_COMM_WORLD);
			MPI_Send(orderMC,N,MPI_INT,0,7,MPI_COMM_WORLD);
			MPI_Send(orderSA,N,MPI_INT,0,8,MPI_COMM_WORLD);
		}
    }



	free(order);
	free(orderSA);
	free(orderMC);
	//free(D);
	free(y);
	free(x);

    MPI_Finalize();


	return 0;
}
