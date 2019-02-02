#include <mpi.h>

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	int p, rank;
	MPI_Comm_size(comm, &p);
	MPI_Comm_rank(comm, &rank);
	
	printf("Hello from rank %i/%i. \n", rank, p);
	
	MPI_Finalize();
	return 0;
}
