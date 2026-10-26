// Emanuele Alterio

#include <mpi.h>
#include <stdio.h>

int somma(int *A, int N)
{

    int NP, myid, token;
    int S, Stot, k;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);

    S = 0;
    for (k = 0; k < N; k++)
    {
        S += A[k];
    }

    Stot = S;
    token = S;

    for (k = 1; k < NP; k++)
    {

        MPI_Send(&token, 1, MPI_INT, (myid + 1) % NP, 10, MPI_COMM_WORLD);
        MPI_Recv(&token, 1, MPI_INT, (myid + NP - 1) % NP, 10, MPI_COMM_WORLD, &status);

        Stot += token;
    }

    return Stot;
}
