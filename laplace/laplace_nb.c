//Emanuele Alterio

#include <stdio.h>
#include <mpi.h>

void laplace(float *A, float *B, float *daprev, float *danext, int N, int LD, int Niter){
    int myid, NP, prev, next;
    int rows, i, j, iter;
    MPI_Request rreq[2], sreq[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &NP);

    rows = N / NP;
    prev = (myid==0) ? MPI_PROC_NULL : myid-1;
    next = (myid == NP-1) ? MPI_PROC_NULL : myid+1;
    
    for(iter=0; iter<Niter; iter++){
        
        MPI_Irecv(daprev, N, MPI_FLOAT, prev, 20, MPI_COMM_WORLD, &rreq[0]);
        MPI_Irecv(danext, N, MPI_FLOAT, next, 10, MPI_COMM_WORLD, &rreq[1]);
        MPI_Isend(&A[0], N, MPI_FLOAT, prev, 10, MPI_COMM_WORLD, &sreq[0]);
        MPI_Isend(&A[(rows-1)*LD], N, MPI_FLOAT, next, 20, MPI_COMM_WORLD, &sreq[1]);

        for(i=1; i<rows-1; i++){
            for(j=1; j<N-1; j++){
                B[i*LD+j] = ( A[(i - 1) * LD + j] + A[(i + 1) * LD + j] + A[i * LD + (j - 1)] + A[i * LD + (j + 1)]) * 0.25f;
            }
        }

        MPI_Wait(&rreq[0], MPI_STATUS_IGNORE);
        if (myid != 0){
            for (j=1; j<N-1; j++){
                B[0 * LD + j] = (daprev[j] + A[1 * LD + j] + A[j-1] + A[j+1]) * 0.25f;
            }
        }

        MPI_Wait(&rreq[1], MPI_STATUS_IGNORE);
        if (myid != NP-1){
            for (j=1; j<N-1; j++){
                B[(rows-1) * LD + j] = (A[(rows-2) * LD + j] + danext[j] + A[(rows-1) * LD + j - 1] + A[(rows-1) * LD + j + 1]) * 0.25f;
            }
        }

        MPI_Wait(&sreq[0], MPI_STATUS_IGNORE);
        MPI_Wait(&sreq[1], MPI_STATUS_IGNORE);

        for (i=1; i<rows-1; i++)
            for (j=1; j<N-1; j++)
                A[i * LD + j] = B[i * LD + j];
 
        if (myid != 0)
            for (j=1; j<N-1; j++)
                A[j] = B[j];
 
        if (myid != NP-1)
            for (j=1; j<N-1; j++)
                A[(rows-1) * LD + j] = B[(rows-1) * LD + j];
    }

}

