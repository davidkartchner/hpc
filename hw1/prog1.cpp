#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>

# define PI 3.14159265358979323846 

int dboard(int N){
// Compute number of darts to be thrown locally
if (rank < N % p) {
    n = (N/p) + 1;
}
else {
    n = N/p;
}

// Throw darts (Generate random coords (x,y))
int m = 0;
for (i=0; i<n; i++){
	// Generate r^2 and theta.  Use a circle with radius 2 for simplicity.
	a = sqrt(2.0f * rand() / (1.0f * RAND_MAX));
    theta = 2.0f*PI*rand()/RAND_MAX;

    // Check whether dart lands in square
    if (fabs(a * cos(theta)) <= 1 && fabs(a * sin(theta)) <= 1){
        m += 1;
    }
}
return m;
}

int main(int argc, char *argv[]){
    // Initialize MPI 
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int p, rank;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &rank);

    // Get passed values for N, R
    if (rank==0) int N = argv[2];
    int R = argv[3];

    // Broadcast N from master processor
    MPI_Bcast(N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int Ms[R] = {};
    for (j=0; j<R; j++){
        MPI_Reduce(&m, &Ms[j], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    // Seed random number generator
    TODO

    // Compute and average pi
    if (rank==0):
        double pi_approx = 0;
        for (j=0; j<R; j++){
             pi_approx +=] = (1.0 * Ms[j])/(1.0 * N);
        }
        pi_approx /= (1.0*R);
        std::cout << pi_approx << std::endl;
	MPI_Finalize();
	return 0;
}