#include <omp.h>

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j, k;
    for (i = 0; i < N1; i++)
        for (j = 0; j < N3; j++)
            for (k = 0; k < N2; k++)
                C[i*ldC + j] = C[i*ldC + j] + A[i*ldA + k] * B[k*ldB + j];
}

void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j, k;
    for (j = 0; j < N3; j++)
        for (i = 0; i < N1; i++)
            for (k = 0; k < N2; k++)
                C[i*ldC + j] = C[i*ldC + j] + A[i*ldA + k] * B[k*ldB + j];
}

void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j, k;
    for (i = 0; i < N1; i++)
        for (k = 0; k < N2; k++)
            for (j = 0; j < N3; j++)
                C[i*ldC + j] = C[i*ldC + j] + A[i*ldA + k] * B[k*ldB + j];
}

void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j, k;
    for (k = 0; k < N2; k++)
        for (i = 0; i < N1; i++)
            for (j = 0; j < N3; j++)
                C[i*ldC + j] = C[i*ldC + j] + A[i*ldA + k] * B[k*ldB + j];
}

void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j, k;
    for (j = 0; j < N3; j++)
        for (k = 0; k < N2; k++)
            for (i = 0; i < N1; i++)
                C[i*ldC + j] = C[i*ldC + j] + A[i*ldA + k] * B[k*ldB + j];
}

void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j, k;
    for (k = 0; k < N2; k++)
        for (j = 0; j < N3; j++)
            for (i = 0; i < N1; i++)
                C[i*ldC + j] = C[i*ldC + j] + A[i*ldA + k] * B[k*ldB + j];
}

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, 
    int N1, int N2, int N3, int dbA, int dbB, int dbC){

        int ii, jj, kk;
        for(ii=0; ii<N1; ii+=dbA){
            for(jj=0; jj<N3; jj+=dbC){
                for(kk=0; kk<N2; kk+=dbB){
                    matmatikj(ldA, ldB, ldC, &A[ii * ldA + kk], &B[kk * ldB + jj], &C[ii * ldC + jj], dbA, dbB, dbC);
                }
            }
        }
}

void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW, int NTCOL){

    int NT = NTROW * NTCOL;

    #pragma omp parallel num_threads(NT)
    {
        int t   = omp_get_thread_num();
        int IDi = t / NTCOL;
        int IDj = t % NTCOL;

        int nr = N1 / NTROW;
        int nc = N3 / NTCOL;

        int r0 = IDi * nr;
        int c0 = IDj * nc;

        matmatblock(ldA, ldB, ldC,
                    &A[r0*ldA],
                    &B[c0],
                    &C[r0*ldC + c0],
                    nr, N2, nc,
                    dbA, dbB, dbC);
    }
}