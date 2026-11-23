#include <stdio.h>
#include <stdlib.h>

#define NMAX 1536
#define LD NMAX
#define DB 256

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);

double get_cur_time();


void set_zero(double *C, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            C[i*LD + j] = 0.0;
}

int main(){
    double *A, *B , *C; 
    int i, j, N;
    double g_ijk, g_jik, g_ikj, g_kij, g_jki, g_kji, g_block;
    double t0, nflops, best_time;
    const char *best_vers;

    A = (double*)malloc((sizeof(double)) * LD * LD);
    B = (double*)malloc((sizeof(double)) * LD * LD);
    C = (double*)malloc((sizeof(double)) * LD * LD);

    for(i=0; i<NMAX; i++){
        for(j=0; j<NMAX; j++){
            A[i*LD+j] = (double)rand() / RAND_MAX;
            B[i*LD+j] = (double)rand() / RAND_MAX;
        }
    }

    printf("Prestazioni in Gflops\n");
    printf("%6s%9s%9s%9s%9s%9s%9s%10s\n", "N", "ijk", "jik", "ikj", "kij", "jki", "kji", "best_vers");

    for (N=256; N<=NMAX; N+=256){
        nflops = 2.0 * N * N * N;

        /* versione ijk */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatijk(LD, LD, LD, A, B, C, N, N, N);
        g_ijk = nflops / (get_cur_time() - t0) / 1.0e9;

        /* versione jik */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatjik(LD, LD, LD, A, B, C, N, N, N);
        g_jik = nflops / (get_cur_time() - t0) / 1.0e9;

        /* versione ikj */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatikj(LD, LD, LD, A, B, C, N, N, N);
        g_ikj = nflops / (get_cur_time() - t0) / 1.0e9;

        /* versione kij */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatkij(LD, LD, LD, A, B, C, N, N, N);
        g_kij = nflops / (get_cur_time() - t0) / 1.0e9;

        /* versione jki */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatjki(LD, LD, LD, A, B, C, N, N, N);
        g_jki = nflops / (get_cur_time() - t0) / 1.0e9;

        /* versione kji */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatkji(LD, LD, LD, A, B, C, N, N, N);
        g_kji = nflops / (get_cur_time() - t0) / 1.0e9;

        best_time = g_ijk;
        best_vers = "ijk";


        if (g_jik > best_time) { best_time = g_jik; best_vers = "jik"; }
        if (g_ikj > best_time) { best_time = g_ikj; best_vers = "ikj"; }
        if (g_kij > best_time) { best_time = g_kij; best_vers = "kij"; }
        if (g_jki > best_time) { best_time = g_jki; best_vers = "jki"; }
        if (g_kji > best_time) { best_time = g_kji; best_vers = "kji"; }

        printf("%6d%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%10s\n", N, g_ijk, g_jik, g_ikj, g_kij, g_jki, g_kji, best_vers);

    }

    /* confronto tra matmatikj e matmatblock */
    printf("\nConfronto ikj / block in Gflops (blocchi %d x %d)\n", DB, DB);
    printf("%6s%9s%9s\n", "N", "ikj", "block");

    for (N=256; N<=NMAX; N+=256){
        nflops = 2.0 * N * N * N;

        /* versione ikj */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatikj(LD, LD, LD, A, B, C, N, N, N);
        g_ikj = nflops / (get_cur_time() - t0) / 1.0e9;

        /* versione a blocchi */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatblock(LD, LD, LD, A, B, C, N, N, N, DB, DB, DB);
        g_block = nflops / (get_cur_time() - t0) / 1.0e9;

        printf("%6d%9.3f%9.3f\n", N, g_ikj, g_block);
    }
    
    free(A); free(B); free(C);
    return 0;
}
