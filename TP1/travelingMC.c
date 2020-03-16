#include <stdlib.h>

void shuffle(int *array, int N) {
    for (int i = 0; i < N - 1; i++) {
        int j = i + rand() / (RAND_MAX / (N - i) + 1);
        
        int t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

void randperm(int *array, int N) {
    for (int i = 0; i < N; i++) {
        array[i] = i;
    }
    
    shuffle(array, N);
}

double MC(double* D, int N, int **town) {
    *town = malloc(sizeof(int) * N);
    double Tdist;
    int i = 0, previous, next1, next2;
    
    randperm(*town, N); // a random permutation of the first n integers
    
    Tdist = D[(*town)[N-1]*N +  (*town)[0]];
    for (int i = 0; i < N-1; i++) { // initial length of route
        Tdist += D[(*town)[i]*N +  (*town)[i + 1]];
    }
    
    while (i < 100) { // stop if no changes occur for 100 trials
        int c = rand() % N; // randomly chooses a town (at position c in route)
        if (c == 0) {
            previous = N - 1;
            next1 = 1;
            next2 = 2;
        } else if (c == N - 2) {
            previous = c - 1;
            next1 = N - 1;
            next2 = 0;
        } else if (c == N - 1) {
            previous = c - 1;
            next1 = 0;
            next2 = 1;
        } else {
            previous = c - 1;
            next1 = c + 1;
            next2 = c + 2;
        }

        // delta=increment in length of route
        double delta =
            D[(*town)[previous]*N +  (*town)[next1]]
                + D[(*town)[c]*N +  (*town)[next2]]
                - D[(*town)[previous]*N +  (*town)[c]]
                - D[(*town)[next1]*N + (*town)[next2]];
        
        // accept or discard change to route
                
        if (delta < 0) {
            // swap order of town(c) and town(c+1) in route  
            int temp = (*town)[c];
            (*town)[c] = (*town)[next1];
            (*town)[next1] = temp;
            
            Tdist += delta;
            i = 0;
        }
        else {
            i++;
        }
    }

    return Tdist;
}

