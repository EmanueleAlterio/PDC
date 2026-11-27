/*
 * Compilazione:   gcc -O3 -fopenmp mainmatmat.c matmat.c c_timer.c -o matmat
 */

#include <stdio.h>
#include <stdlib.h>

#define NMAX 1536
#define NMAX_THREAD 2048
#define LD NMAX_THREAD
#define DB 256

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW, int NTCOL);

double get_cur_time();

void set_zero(double *C, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            C[i * LD + j] = 0.0;
}

/* configurazioni della griglia di thread richieste dalla slide */
int ntrow[4] = {1, 1, 2, 2};
int ntcol[4] = {1, 2, 2, 4};

/* test di correttezza di matmatthread: N = 8, blocchi 2x2, confronto con matmatikj */
void test_correttezza()
{
    double A[64], B[64], C[64], R[64];
    double diff, d;
    int i, c;

    for (i = 0; i < 64; i++)
    {
        A[i] = (double)rand() / RAND_MAX;
        B[i] = (double)rand() / RAND_MAX;
        R[i] = 0.0;
    }

    /* risultato di riferimento */
    matmatikj(8, 8, 8, A, B, R, 8, 8, 8);

    printf("Test di correttezza matmatthread (N = 8, blocchi 2 x 2)\n");
    for (c = 0; c < 4; c++)
    {
        for (i = 0; i < 64; i++)
            C[i] = 0.0;
        matmatthread(8, 8, 8, A, B, C, 8, 8, 8, 2, 2, 2, ntrow[c], ntcol[c]);

        diff = 0.0;
        for (i = 0; i < 64; i++)
        {
            d = C[i] - R[i];
            if (d < 0)
                d = -d;
            if (d > diff)
                diff = d;
        }
        printf("  griglia %d x %d : diff max = %.2e  %s\n",
               ntrow[c], ntcol[c], diff, diff < 1e-12 ? "OK" : "ERRORE");
    }
}

int main()
{
    double *A, *B, *C;
    int i, j, N;
    double g_ijk, g_jik, g_ikj, g_kij, g_jki, g_kji, g_block;
    double t0, nflops, best_time;
    double t_block, t_thread, g_thread, speedup, efficienza;
    int c, NT;
    const char *best_vers;

    A = (double *)malloc((sizeof(double)) * LD * LD);
    B = (double *)malloc((sizeof(double)) * LD * LD);
    C = (double *)malloc((sizeof(double)) * LD * LD);

    for (i = 0; i < LD; i++)
    {
        for (j = 0; j < LD; j++)
        {
            A[i * LD + j] = (double)rand() / RAND_MAX;
            B[i * LD + j] = (double)rand() / RAND_MAX;
        }
    }

    /*
    printf("Prestazioni in Gflops\n");
    printf("%6s%9s%9s%9s%9s%9s%9s%10s\n", "N", "ijk", "jik", "ikj", "kij", "jki", "kji", "best_vers");

    for (N = 256; N <= NMAX; N += 256)
    {
        nflops = 2.0 * N * N * N;

        // versione ijk
        set_zero(C, N);
        t0 = get_cur_time();
        matmatijk(LD, LD, LD, A, B, C, N, N, N);
        g_ijk = nflops / (get_cur_time() - t0) / 1.0e9;

        // versione jik
        set_zero(C, N);
        t0 = get_cur_time();
        matmatjik(LD, LD, LD, A, B, C, N, N, N);
        g_jik = nflops / (get_cur_time() - t0) / 1.0e9;

        // versione ikj
        set_zero(C, N);
        t0 = get_cur_time();
        matmatikj(LD, LD, LD, A, B, C, N, N, N);
        g_ikj = nflops / (get_cur_time() - t0) / 1.0e9;

        // versione kij
        set_zero(C, N);
        t0 = get_cur_time();
        matmatkij(LD, LD, LD, A, B, C, N, N, N);
        g_kij = nflops / (get_cur_time() - t0) / 1.0e9;

        // versione jki
        set_zero(C, N);
        t0 = get_cur_time();
        matmatjki(LD, LD, LD, A, B, C, N, N, N);
        g_jki = nflops / (get_cur_time() - t0) / 1.0e9;

        // versione kji
        set_zero(C, N);
        t0 = get_cur_time();
        matmatkji(LD, LD, LD, A, B, C, N, N, N);
        g_kji = nflops / (get_cur_time() - t0) / 1.0e9;

        best_time = g_ijk;
        best_vers = "ijk";

        if (g_jik > best_time)
        {
            best_time = g_jik;
            best_vers = "jik";
        }
        if (g_ikj > best_time)
        {
            best_time = g_ikj;
            best_vers = "ikj";
        }
        if (g_kij > best_time)
        {
            best_time = g_kij;
            best_vers = "kij";
        }
        if (g_jki > best_time)
        {
            best_time = g_jki;
            best_vers = "jki";
        }
        if (g_kji > best_time)
        {
            best_time = g_kji;
            best_vers = "kji";
        }

        printf("%6d%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%10s\n", N, g_ijk, g_jik, g_ikj, g_kij, g_jki, g_kji, best_vers);
    }

    // confronto tra matmatikj e matmatblock
    printf("\nConfronto ikj / block in Gflops (blocchi %d x %d)\n", DB, DB);
    printf("%6s%9s%9s\n", "N", "ikj", "block");

    for (N = 256; N <= NMAX; N += 256)
    {
        nflops = 2.0 * N * N * N;

        // versione ikj
        set_zero(C, N);
        t0 = get_cur_time();
        matmatikj(LD, LD, LD, A, B, C, N, N, N);
        g_ikj = nflops / (get_cur_time() - t0) / 1.0e9;

        // versione a blocchi
        set_zero(C, N);
        t0 = get_cur_time();
        matmatblock(LD, LD, LD, A, B, C, N, N, N, DB, DB, DB);
        g_block = nflops / (get_cur_time() - t0) / 1.0e9;

        printf("%6d%9.3f%9.3f\n", N, g_ikj, g_block);
    }
    */

    /* ------------------- matmatthread ------------------- */
    printf("\n");
    test_correttezza();

    /* confronto tra matmatblock e matmatthread per N = 1024, 2048 */
    for (N = 1024; N <= NMAX_THREAD; N *= 2)
    {
        nflops = 2.0 * N * N * N;

        /* riferimento: matmatblock (sequenziale) */
        set_zero(C, N);
        t0 = get_cur_time();
        matmatblock(LD, LD, LD, A, B, C, N, N, N, DB, DB, DB);
        t_block = get_cur_time() - t0;
        g_block = nflops / t_block / 1.0e9;

        printf("\nConfronto block / thread, N = %d (blocchi %d x %d)\n", N, DB, DB);
        printf("%-12s%6s%10s%10s%10s%12s\n", "versione", "NT", "tempo(s)", "Gflops", "SpeedUp", "Efficienza");
        printf("%-12s%6d%10.3f%10.3f%10s%12s\n", "matmatblock", 1, t_block, g_block, "-", "-");

        /* matmatthread con le 4 griglie di thread */
        for (c = 0; c < 4; c++)
        {
            NT = ntrow[c] * ntcol[c];

            set_zero(C, N);
            t0 = get_cur_time();
            matmatthread(LD, LD, LD, A, B, C, N, N, N, DB, DB, DB, ntrow[c], ntcol[c]);
            t_thread = get_cur_time() - t0;

            g_thread = nflops / t_thread / 1.0e9;
            speedup = t_block / t_thread; /* rispetto a matmatblock */
            efficienza = speedup / NT;

            printf("thread %dx%d %6d%10.3f%10.3f%10.2f%12.2f\n",
                   ntrow[c], ntcol[c], NT, t_thread, g_thread, speedup, efficienza);
        }
    }

    free(A);
    free(B);
    free(C);
    return 0;
}