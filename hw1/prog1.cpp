#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>

# define PI 3.14159265358979323846 

int dboard(int N, int rank, int p){
    // Compute number of darts to be thrown locally
    int n = N/p;
    if (rank < N % p) {
        n += 1;
    }

    // Throw darts (Generate random coords (x,y))
    int m = 0;
    for (int i=0; i<n; i++){
        // Generate r^2 and theta.  Use a circle with radius 2 for simplicity.
        float a = sqrt(2.0f * rand() / (1.0f * RAND_MAX));
        float theta = 2.0f*PI*rand()/RAND_MAX;

        // Check whether dart lands in square
        if (fabs(a * cos(theta)) <= 1 && fabs(a * sin(theta)) <= 1){
            m += 1;
        }
    }
    return m;
}

int main(int argc, char *argv[]){
    //Start timing program
    float tic = time(NULL);

    // Initialize MPI 
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int p, rank;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &rank);

    // Get passed values for N, R
    int N = 0;
    if (rank==0) N = *argv[-2];
    int R = *argv[-1];

    // Broadcast N from master processor
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int Ms[R];
    for (int j=0; j<R; j++){
        // Seed random number generator
        srand(rank);
        int m = dboard(N, rank, p);
        MPI_Reduce(&m, &Ms[j], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Compute and average pi
    double pi_approx = 0;
    if (rank==0){
        for (int j=0; j<R; j++){
             pi_approx += (1.0 * Ms[j])/(1.0 * N);
        }
        pi_approx /= (1.0*R);
        std::cout << pi_approx << std::endl;
    }
	MPI_Finalize();

    // Print total compute time
    float toc = time(NULL);
    float seconds = difftime(toc, tic);
    std::cout << "Total Computation time: " << seconds << " Seconds" << std::endl;

	return 0;
}