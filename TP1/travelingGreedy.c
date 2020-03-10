
#include <math.h>


double greedy(double* D, int N, int** order, int starting_town){

	double distance = 0;
	*order = (int*)malloc(sizeof(int)*N);
	int* visited = malloc(sizeof(int)*N);

	// Visited towns will be recorded with a value of 1 in its corresponing index of the visited array
	for (int i = 0; i < N; ++i)
	{
		visited[i] = 0;
	}


	(*order)[0] = starting_town;
	visited[starting_town] = 1;

	for (int i = 0; i < N-1; i++){
		double min = INFINITY;
		for (int j = 0; j < N; j++){
			
			// Check if town has not yet been visited
			if (!visited[j]){
				// Check if town is currently the closest
				if(D[(*order)[i]*N + j] < min){
					min = D[(*order)[i]*N + j];
					(*order)[i+1] = j;
				}
			} 
		}

		visited[(*order)[i+1]] = 1;
		distance += D[(*order)[i]*N + (*order)[i+1]];
	}

	distance += D[(*order)[N-1]*N + (*order)[0]];

	free(visited);

	return distance;

}