#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

extern double get_cur_time();

double maxsum(int N, int ld, double *matrix, int NT);

double** matrix_allocation(size_t N){
    size_t i;

    double *pool = malloc(N * N * sizeof(double));
    if (pool == NULL){
        printf("Error in matrix allocation.\n");
        return NULL;
    }

    double **matrix = malloc(N * sizeof(double*));
    if (matrix == NULL){
        printf("Error in matrix allocation.\n");
        free(pool);
        return NULL;
    }

    for(i=0; i<N; i++){
        matrix[i] = &pool[i*N];
    }

    for(i=0; i<N*N; i++){
        pool[i] = rand() % 100;
    }

    printf("%zux%zu Matrix has been allocated.\n", N, N);

    return matrix;
}

void matrix_free(double** matrix){
    if (matrix==NULL){
        printf("Nothing to free.\n");
        return;
    }

    free(matrix[0]);
    free(matrix);
    printf("Memory free correctly.\n");
}

int main(void){
    int ld = 5400;

    double **matrix = matrix_allocation((size_t)ld);
    if (matrix == NULL){
        return 1;
    }

    int size_mat[] = {800, 1600, 2400, 3200, 4000, 4800, 5400};
    int count_size = sizeof(size_mat) / sizeof(size_mat[0]);

    int NT[] = {1, 2, 4, 8};
    int count_NT = sizeof(NT) / sizeof(NT[0]);

    for(int i=0; i<count_size; i++){
        printf("MaxSum on size %dx%d;\n", size_mat[i], size_mat[i]);
        double t_seq = 0.0;

        for(int j=0; j<count_NT; j++){
            double start = get_cur_time();
            double max = maxsum(size_mat[i], ld, matrix[0], NT[j]);
            double end = get_cur_time();
            double time_execution = end - start;

            if (NT[j] == 1){
                t_seq = time_execution;
            }

            double speedup = t_seq / time_execution;
            double efficiency = speedup / NT[j];

            printf("NT=%d | max=%f | time=%.6fs | speedup=%.3f | efficiency=%.3f\n",
                   NT[j], max, time_execution, speedup, efficiency);
        }
    }

    matrix_free(matrix);

    return 0;
}