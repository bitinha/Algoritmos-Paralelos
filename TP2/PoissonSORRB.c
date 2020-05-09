
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


void plotGraph(float tol, int N, float *w) {
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    int i;
    
    fprintf(gnuplotPipe, "unset key\n");
    fprintf(gnuplotPipe, "set style line 11 lc rgb '#808080' lt 1\n");
    fprintf(gnuplotPipe, "set border 3 front ls 11\n");
    fprintf(gnuplotPipe, "set tics nomirror out scale 0.75\n");
    fprintf(gnuplotPipe, "set cbtics scale 0\n");
    fprintf(gnuplotPipe, "load 'parula.pal'\n");
    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set yrange[0:%d]\n",N);
    fprintf(gnuplotPipe, "set xrange[0:%d]\n",N);
    fprintf(gnuplotPipe, "set output '%d,%f.png'\n", N, tol);
    fprintf(gnuplotPipe, "plot '-' u ($1):($2):($3) matrix with image\n");

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			fprintf(gnuplotPipe, "%f\t", w[i*N+j]);

		}
		fprintf(gnuplotPipe,"\n");
	}

    fflush(gnuplotPipe);
}




float maxDifference(float *a, float *b, int N){
	float max = fabsf(a[0]-b[0]);
	float m;
	for (int i = 1; i < N; ++i)
	{
		m = fabsf(a[i]-b[i]);
		if (m > max){
			max = m;
		}
	}
	return max;
}


int main(int argc, char const *argv[])
{


	if (argc < 2)
	{
		printf("Insira um valor para N seguido de um valor para a tolerancia\n");
		return 1;
	}

	int N = atoi(argv[1]);
	float tol = atof(argv[2]);

	float *w = (float *)malloc(N*N*sizeof(float));


	for (int j = 0; j < N; ++j)
	{
		w[j] = 100.0;
	}
	for (int i = 1; i < N-1; ++i)
	{
		w[i*N] = 100.0;
		for (int j = 1; j < N-1; ++j)
		{
			w[i*N + j] = 50.0;
		}
		w[(i+1)*N-1] = 100.0;
	}
	for (int j = 0; j < N; ++j)
	{
		w[N*(N-1)+j] = 0.0;
	}


	float p = 2/(1+sin(M_PI/(N-1)));

	float diff = tol + 1;

	int iter = 0;
	float *u = (float *)malloc(N*N*sizeof(float));

	while(diff > tol)
	{
		memcpy(u,w, N*N*sizeof(float));
		for (int i = 1; i < N-1; ++i)
		{
			for (int j = 1 + i%2; j < N-1; j+=2)
			{
				w[i*N+j]=(1-p)*w[i*N+j]+p*(w[(i-1)*N+j]+w[i*N+j-1]+w[i*N+j+1]+w[(i+1)*N+j])/4;
			}
		}

		for (int i = 1; i < N-1; ++i)
		{
			for (int j = 1 + (i+1)%2; j < N-1; j+=2)
			{
				w[i*N+j]=(1-p)*w[i*N+j]+p*(w[(i-1)*N+j]+w[i*N+j-1]+w[i*N+j+1]+w[(i+1)*N+j])/4;
			}
		}
		
		iter++;
		diff = maxDifference(w,u,N*N);

		//printf("%f\n", diff);
	}

	plotGraph(tol,N,w);

	// for(int i = 0; i < N; i++){
	// 	for(int j = 0; j < N; j++){
	// 		printf("%f\t", w[i*N+j]);
	// 	}
	// 	printf("\n");
	// }

	return 0;
}