#include <stdio.h>
#include <mpi.h>

#define N 10

int somma(int*,int);

int main (int argc, char *argv[]){

    int myid;
    int A[N], Stot, i;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    for(i=0; i<N; i++){
        A[i] = myid * 10 + i;
    }

    Stot = somma(A,N);

    printf("Processo %d: somma totale = %d\n", myid, Stot);


    MPI_Finalize();

    return 0;

}